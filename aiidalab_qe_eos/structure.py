import ipywidgets as ipw
import traitlets as tl
from ase import Atoms

class EOSStructure(ipw.VBox):
    """Create structure use ASE builder."""

    structure = tl.Instance(Atoms, allow_none=True)

    def __init__(self, title="EOS"):
        self.title = title
        self.symbol = ipw.Text(placeholder="Au")
        self.create_structure_btn = ipw.Button(
            description="Generate bulk",
            button_style="primary",
            tooltip="Generate bulk from symbol",
        )
        self.create_structure_btn.on_click(self._on_button_pressed)

        super().__init__(
            [ipw.HBox([self.symbol, self.create_structure_btn])]
        )

    def _on_button_pressed(self, change=None):
        """Create a ase structure when button is pressed."""
        from ase.build import bulk
        if not self.symbol.value:
            return
        self.structure = bulk(self.symbol.value)


    @tl.default("structure")
    def _default_structure(self):
        return None
