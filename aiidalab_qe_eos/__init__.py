from aiidalab_qe_eos.setting import Setting
from aiidalab_qe_eos.workchain import workchain_and_builder
from aiidalab_qe_eos.result import Result
from aiidalab_qe.common.panel import OutlinePanel


class EosOutline(OutlinePanel):
    title = "Equation of State (EOS)"


eos ={
"outline": EosOutline,
"setting": Setting,
"workchain": workchain_and_builder,
"result": Result,
}
