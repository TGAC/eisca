
from dominate.tags import div, h1, h4, p, code
from dominate.tags import figure, img
import base64
from ezcharts.layout.snippets import DataTable, Grid, Tabs
from ezcharts.layout.resource import (ImageResource, ScriptResource, StyleResource, ThemeResource)
from ezcharts.layout.util import inline, load_json, transpile, cls, css
from ezcharts.layout.base import IClasses, IStyles, Snippet
from ezcharts.layout.snippets.banner import IBannerClasses, IBannerStyles, IBadgeClasses, IBadgeStyles, Badge
from ezcharts.components.reports.labs import ILabsAddendumClasses



class EIBanner(Snippet):
    """A styled div tag containing a heading and badges."""

    TAG = 'div'

    def __init__(
        self,
        report_title: str,
        workflow_name: str,
        styles: IBannerStyles = IBannerStyles(),
        classes: IBannerClasses = IBannerClasses(),
        default_content: bool = True
    ) -> None:
        """Create styled banner."""
        super().__init__(
            styles=styles,
            classes=classes,
            className=classes.container,
            style=styles.container)

        with self:
            with div(className=classes.inner, style=styles.inner):
                if not default_content:
                    return
                h1(report_title)
                p(
                    f"Results generated through the {workflow_name} nextflow "
                    "workflow developed by Earlham Institute.",
                    className="py-3 fs-5")
                self.badges = div(className="d-flex flex-wrap")

    def add_badge(
        self,
        title: str,
        bg_class=None,
        styles: IBadgeStyles = IBadgeStyles(),
        classes: IBadgeClasses = IBadgeClasses(),
    ) -> None:
        """Add a badge to the banner."""
        with self.badges:
            Badge(title, styles=styles, classes=classes, bg_class=bg_class)


class EILabsAddendum(Snippet):
    """A styled footer component for use in a Report."""

    TAG = 'div'

    def __init__(
        self,
        workflow_name: str,
        workflow_version: str,
        classes: ILabsAddendumClasses = ILabsAddendumClasses(),
        use_defaults: bool = True
    ) -> None:
        """Create tag."""
        super().__init__(
            styles=None,
            classes=classes,
            className=classes.container)

        with self:
            self.container = div(className=classes.inner)

        if use_defaults:
            with self.container:
                h4('About this report', className="pb-3")
                p(
                    "This report was produced using the ",
                    code(f"TGAC/{workflow_name}"),
                    f" nextflow workflow ({workflow_version})."
                )
                # p(
                #     "Oxford Nanopore Technologies products are not "
                #     "intended for use for health assessment or to "
                #     "diagnose, treat, mitigate, cure or prevent any "
                #     "disease or condition.")



def plots_from_image_files(path, meta=None, ncol=1, suffix=['.png'], widths=['1200']):
    """Create a plots section which uses image files.
    PARAMS:
        meta: can be sample or group, by which images are shown in tabs
        ncol: the number of columns for Grid
        suffix: list of suffix used to match the image files
        widths: the column withs of Grid columns
    """

    if meta in ['sample', 'group']:
        tabs = Tabs()
        for folder in sorted([f for f in path.iterdir() if f.is_dir()]):
            if folder.name.startswith(f'{meta}_'):
                sample_id = folder.name.split('_', 1)[1]
                with tabs.add_tab(sample_id):
                    with Grid(columns=ncol):
                        pathes = [next(folder.glob('*'+sf)) for sf in suffix]
                        for img_path, width in zip(pathes, widths):
                            with open(img_path, 'rb') as fh:
                                b64img = base64.b64encode(fh.read()).decode()
                                with figure(cls='text-center'):
                                    img(src=f'data:image/png;base64,{b64img}', width=width)
    else:
        with Grid(columns=ncol):
            pathes = [next(path.glob('*'+sf)) for sf in suffix]
            for img_path, width in zip(pathes, widths):
                with open(img_path, 'rb') as fh:
                    b64img = base64.b64encode(fh.read()).decode()
                    with figure(cls='text-center'):
                        img(src=f'data:image/png;base64,{b64img}', width=width)        


class EPI2MELabsLogo(div):
    """Labs logo element."""

    def __init__(self) -> None:
        """Create a div with an SVG logo inside."""
        super().__init__(
            inline(ImageResource('LAB_logo.svg').data_file),
            tagname='div',
            style="width: 35px; height: 35px;",
            className="d-flex",
            alt="EPI2ME Labs Logo")