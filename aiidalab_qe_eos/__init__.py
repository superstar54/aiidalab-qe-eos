from aiidalab_qe_eos.setting import Setting
from aiidalab_qe_eos.workchain import workchain_and_builder
from aiidalab_qe_eos.result import Result
from aiidalab_qe.panel import OutlinePanel


class Outline(OutlinePanel):
    title = "Equation of State (EOS)"


property ={
"outline": Outline,
"setting": Setting,
"workchain": workchain_and_builder,
"result": Result,
}
