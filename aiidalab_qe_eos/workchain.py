"""Implementation of the MultiplyAddWorkChain for testing and demonstration purposes."""
from aiida.common import AttributeDict
from aiida.engine import ToContext, WorkChain, calcfunction
from aiida.orm import AbstractCode, Int, Float, Dict, Code, StructureData, load_code
from aiida.plugins import WorkflowFactory
from aiida_quantumespresso.utils.mapping import prepare_process_inputs


PwBaseWorkChain = WorkflowFactory('quantumespresso.pw.base')


class EOSWorkChain(WorkChain):
    """WorkChain to calcalculate the equation of state of a crystal."""

    @classmethod
    def define(cls, spec):
        """Specify inputs and outputs."""
        super().define(spec)
        spec.expose_inputs(
            PwBaseWorkChain,
            namespace='scf',
            exclude=('clean_workdir', 'pw.structure', 'pw.parent_folder'),
            namespace_options={
                'required': False,
                'populate_defaults': False,
                'help': 'Inputs for the `PwBaseWorkChain` of the `scf` calculation.',
            }
        )
        spec.input('structure', valid_type=StructureData)
        spec.input('scale', valid_type=Float)
        spec.input('npoint', valid_type=Int)
        spec.outline(
            cls.run_eos,
            cls.result,
        )
        spec.output('eos', valid_type=Dict)
        spec.exit_code(400, 'ERROR_NEGATIVE_NUMBER', message='The result is a negative number.')

    @classmethod
    def get_builder_from_protocol(
        cls,
        codes=None,
        structure=None,
        parameters=None,
    ):
        builder = cls.get_builder()
        builder.structure = structure
        protocol = parameters["basic"].get('protocol', "fast")
        scf_parameters, scf_overrides = cls.get_scf_parameters(parameters)
        # scf
        pw_code = load_code(codes.get('pw_code'))
        args = (pw_code, structure, protocol)
        scf = PwBaseWorkChain.get_builder_from_protocol(
            *args, 
            overrides=scf_overrides,
            **scf_parameters,
        )
        builder.scf = scf
        builder.scale = parameters["eos"]["scale"]
        builder.npoint = parameters["eos"]["npoint"]
        return builder

    @classmethod
    def get_scf_parameters(cls, parameters):
        # developer should get the plugin parameters and override from the parameters
        new_parameters = parameters["basic"]
        pw = parameters["advance"].get("pw", {})
        pw["pseudo_family"] = parameters["advance"].get("pseudo_family", None)
        overrides = {
            "pw": pw,
        }

        return new_parameters, overrides

    def run_eos(self):
        """Run all scf calculations."""
        import numpy as np
        factors = {f"s{x}":x for x in np.linspace(1-self.inputs.scale.value, 1+self.inputs.scale.value, self.inputs.npoint.value)}
        self.ctx.factors = factors
        structure = self.inputs.structure
        futures = {}
        for label, factor in factors.items():
            atoms = structure.get_ase()
            atoms.set_cell(atoms.get_cell()*factor, scale_atoms = True)
            scaled_structure = StructureData(ase=atoms)
            inputs = AttributeDict(self.exposed_inputs(PwBaseWorkChain, namespace="scf"))
            inputs.metadata.call_link_label = "pw"
            inputs.pw.structure = scaled_structure
            inputs.pw.parameters = inputs.pw.parameters.get_dict()
            inputs = prepare_process_inputs(PwBaseWorkChain, inputs)
            running = self.submit(PwBaseWorkChain, **inputs)
            futures[label] = running
            self.report(f"Running an SCF calculation with scale factor {factor}")
        self.report(f"jobs: {futures}")
        return ToContext(**futures)


    def result(self):
        """Add the result to the outputs."""
        vs = []
        es = []
        self.report("factors: ", self.ctx.factors)
        for label, factor in self.ctx.factors.items():
            result = self.ctx[label].outputs.output_parameters
            vs.append(result.dict.volume)
            es.append(result.dict.energy)
            unit = result.dict.energy_units
        eos = {"volume": vs, "energy": es, "unit": unit}
        self.out("eos", Dict(eos))