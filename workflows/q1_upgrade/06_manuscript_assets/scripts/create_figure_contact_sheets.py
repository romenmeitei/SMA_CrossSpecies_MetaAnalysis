from pathlib import Path
import sys

from PIL import Image, ImageDraw, ImageFont

root = Path(sys.argv[1]).resolve()
audit = root / "audit"

audit.mkdir(
    parents=True,
    exist_ok=True,
)


def make_sheet(paths, output, title, columns=2):
    paths = list(paths)

    if not paths:
        raise RuntimeError(
            f"No PNG figures found for {title}"
        )

    cell_width = 1400
    cell_height = 950
    label_height = 75
    margin = 45
    title_height = 110

    rows = (
        len(paths) + columns - 1
    ) // columns

    canvas = Image.new(
        "RGB",
        (
            columns * cell_width
            + (columns + 1) * margin,
            title_height
            + rows * (
                cell_height + label_height
            )
            + (rows + 1) * margin,
        ),
        "white",
    )

    draw = ImageDraw.Draw(canvas)

    try:
        title_font = ImageFont.truetype(
            "DejaVuSans-Bold.ttf",
            34,
        )
        label_font = ImageFont.truetype(
            "DejaVuSans.ttf",
            22,
        )
    except OSError:
        title_font = ImageFont.load_default()
        label_font = ImageFont.load_default()

    draw.text(
        (margin, 35),
        title,
        fill="black",
        font=title_font,
    )

    for index, path in enumerate(paths):
        row = index // columns
        column = index % columns

        x0 = margin + column * (
            cell_width + margin
        )

        y0 = title_height + margin + row * (
            cell_height
            + label_height
            + margin
        )

        with Image.open(path) as image:
            image = image.convert("RGB")
            image.thumbnail(
                (
                    cell_width,
                    cell_height,
                ),
                Image.Resampling.LANCZOS,
            )

            x = x0 + (
                cell_width - image.width
            ) // 2

            y = y0 + (
                cell_height - image.height
            ) // 2

            canvas.paste(image, (x, y))

        draw.rectangle(
            (
                x0,
                y0,
                x0 + cell_width,
                y0 + cell_height,
            ),
            outline="#B0B0B0",
            width=2,
        )

        draw.text(
            (
                x0,
                y0 + cell_height + 15,
            ),
            path.name,
            fill="black",
            font=label_font,
        )

    canvas.save(
        output,
        dpi=(150, 150),
    )

    print("Created:", output)


make_sheet(
    sorted(
        (root / "figures" / "main")
        .glob("*.png")
    ),
    audit / "main_figures_contact_sheet.png",
    "SMA manuscript main figures",
)

make_sheet(
    sorted(
        (root / "figures" / "supplementary")
        .glob("*.png")
    ),
    audit
    / "supplementary_figures_contact_sheet.png",
    "SMA manuscript supplementary figures",
)
