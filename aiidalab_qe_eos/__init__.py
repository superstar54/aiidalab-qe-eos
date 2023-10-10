from .setting import Setting
from .workchain import workchain_and_builder
from .result import Result
from .structure import EOSStructure
from aiidalab_qe.common.panel import OutlinePanel


class EosOutline(OutlinePanel):
    title = "Equation of State (EOS)"


eos ={
"outline": EosOutline,
"importer": EOSStructure,
"setting": Setting,
"workchain": workchain_and_builder,
"result": Result,
}
