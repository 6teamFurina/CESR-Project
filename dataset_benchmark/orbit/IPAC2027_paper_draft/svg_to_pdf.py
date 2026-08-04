#!/usr/bin/env python3
"""Convert the three simple project SVG figures to PDF for LaTeX inclusion.

The maintained figures use only rect, line, polyline, polygon, circle, and text
elements.  Keeping this converter local avoids changing the scientific plotting
code or requiring an external SVG application in the paper-build environment.
"""

from __future__ import annotations

import re
import xml.etree.ElementTree as ET
from pathlib import Path

from reportlab.graphics import renderPDF
from reportlab.graphics.shapes import Circle, Drawing, Group, Line, Polygon, PolyLine, Rect, String
from reportlab.lib import colors
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont


HERE = Path(__file__).resolve().parent
FIGURES = HERE / "figures"
SOURCES = (
    "vertical_parity_crossover_loglog.svg",
    "vertical_parity_growth_linear.svg",
    "vertical_parity_odd_fraction_percentiles.svg",
)


def register_fonts() -> tuple[str, str]:
    candidates = (
        (
            Path(r"C:\Windows\Fonts\arial.ttf"),
            Path(r"C:\Windows\Fonts\arialbd.ttf"),
        ),
        (
            Path(r"C:\Windows\Fonts\calibri.ttf"),
            Path(r"C:\Windows\Fonts\calibrib.ttf"),
        ),
    )
    for regular, bold in candidates:
        if regular.exists() and bold.exists():
            pdfmetrics.registerFont(TTFont("ProjectSans", regular))
            pdfmetrics.registerFont(TTFont("ProjectSans-Bold", bold))
            return "ProjectSans", "ProjectSans-Bold"
    return "Helvetica", "Helvetica-Bold"


FONT_REGULAR, FONT_BOLD = register_fonts()


def tag_name(element: ET.Element) -> str:
    return element.tag.rsplit("}", 1)[-1]


def number(value: str | None, total: float | None = None) -> float:
    if value is None:
        return 0.0
    value = value.strip()
    if value.endswith("%"):
        if total is None:
            raise ValueError(f"Percentage length lacks a reference: {value}")
        return total * float(value[:-1]) / 100.0
    match = re.match(r"[-+0-9.eE]+", value)
    if not match:
        raise ValueError(f"Unsupported SVG length: {value}")
    return float(match.group(0))


def paint(value: str | None, opacity: float = 1.0):
    if value in (None, "none"):
        return None
    base = colors.toColor(value)
    if opacity >= 1.0:
        return base
    return colors.Color(base.red, base.green, base.blue, alpha=opacity)


def style_value(style: dict[str, str], name: str, default: str | None = None) -> str | None:
    return style.get(name, default)


def apply_stroke(shape, style: dict[str, str]) -> None:
    shape.strokeColor = paint(style_value(style, "stroke", "none"))
    shape.strokeWidth = number(style_value(style, "stroke-width", "1"))
    dash = style_value(style, "stroke-dasharray")
    if dash:
        shape.strokeDashArray = [float(item) for item in re.findall(r"[-+0-9.eE]+", dash)]


def apply_fill(shape, style: dict[str, str]) -> None:
    opacity = float(style_value(style, "fill-opacity", "1"))
    shape.fillColor = paint(style_value(style, "fill", "black"), opacity)


def point_values(value: str, height: float) -> list[float]:
    raw = [float(item) for item in re.findall(r"[-+0-9.eE]+", value)]
    if len(raw) % 2:
        raise ValueError(f"Odd coordinate count in points={value!r}")
    converted: list[float] = []
    for index in range(0, len(raw), 2):
        converted.extend((raw[index], height - raw[index + 1]))
    return converted


def convert_element(
    element: ET.Element,
    drawing: Drawing,
    width: float,
    height: float,
    inherited: dict[str, str],
) -> None:
    style = dict(inherited)
    style.update(element.attrib)
    name = tag_name(element)

    if name in {"svg", "g"}:
        for child in element:
            convert_element(child, drawing, width, height, style)
        return
    if name in {"title", "desc"}:
        return

    if name == "rect":
        x = number(element.get("x"), width)
        y = number(element.get("y"), height)
        rect_width = number(element.get("width"), width)
        rect_height = number(element.get("height"), height)
        shape = Rect(x, height - y - rect_height, rect_width, rect_height)
        apply_fill(shape, style)
        apply_stroke(shape, style)
        drawing.add(shape)
        return

    if name == "line":
        shape = Line(
            number(element.get("x1")),
            height - number(element.get("y1")),
            number(element.get("x2")),
            height - number(element.get("y2")),
        )
        apply_stroke(shape, style)
        drawing.add(shape)
        return

    if name in {"polyline", "polygon"}:
        points = point_values(element.get("points", ""), height)
        shape = PolyLine(points) if name == "polyline" else Polygon(points)
        if name == "polygon":
            apply_fill(shape, style)
        apply_stroke(shape, style)
        drawing.add(shape)
        return

    if name == "circle":
        shape = Circle(
            number(element.get("cx")),
            height - number(element.get("cy")),
            number(element.get("r")),
        )
        apply_fill(shape, style)
        apply_stroke(shape, style)
        drawing.add(shape)
        return

    if name == "text":
        x = number(element.get("x"))
        y_svg = number(element.get("y"))
        anchor = style_value(style, "text-anchor", "start")
        font_size = number(style_value(style, "font-size", "12"))
        weight = style_value(style, "font-weight", "normal")
        font_name = FONT_BOLD if weight in {"bold", "600", "700"} else FONT_REGULAR
        content = "".join(element.itertext())
        text = String(
            x,
            height - y_svg,
            content,
            fontName=font_name,
            fontSize=font_size,
            fillColor=paint(style_value(style, "fill", "black")),
            textAnchor=anchor,
        )
        transform = element.get("transform", "")
        rotation = re.fullmatch(
            r"rotate\(\s*(-?90(?:\.0+)?)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s*\)",
            transform,
        )
        if rotation:
            angle = float(rotation.group(1))
            cx = float(rotation.group(2))
            cy = height - float(rotation.group(3))
            text.x = 0
            text.y = 0
            group = Group()
            group.translate(cx, cy)
            group.rotate(-angle)
            group.add(text)
            drawing.add(group)
        else:
            drawing.add(text)
        return

    raise ValueError(f"Unsupported SVG element <{name}>")


def convert(source: Path) -> Path:
    root = ET.parse(source).getroot()
    width = number(root.get("width"))
    height = number(root.get("height"))
    drawing = Drawing(width, height)
    convert_element(root, drawing, width, height, {})
    destination = source.with_suffix(".pdf")
    renderPDF.drawToFile(drawing, str(destination))
    return destination


def main() -> int:
    for filename in SOURCES:
        destination = convert(FIGURES / filename)
        print(destination)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
