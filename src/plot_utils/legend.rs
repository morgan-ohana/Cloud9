use plotters::coord::Shift;
use plotters::prelude::*;
use plotters::style::text_anchor::{HPos, Pos, VPos};
use std::error::Error;

pub enum LegendMarker<DB: DrawingBackend> {
    Line(RGBColor),
    Dashed(RGBColor, i32, i32), // (color, dash_on, dash_off) in px
    /// (area, y_center, x0, x1) — for exotic entries like your multi-color dashed line
    Custom(Box<dyn Fn(&DrawingArea<DB, Shift>, i32, i32, i32) -> Result<(), Box<dyn Error>>>),
}

pub struct LegendEntry<DB: DrawingBackend> {
    pub tex: String,   // LaTeX source — typeset by the pdf_tex pipeline
    pub proxy: String, // plain-text approximation of the typeset label, sizes the box
    pub marker: LegendMarker<DB>,
}

pub enum LegendAnchor {
    UpperRight,
    UpperLeft,
}

/// Approximate the *typeset* label from its LaTeX source, for width estimation.
/// Only needs to be close — DejaVu (used for measuring) is wider than Latin
/// Modern Sans (used by LaTeX), so errors lean toward a generous box.
pub fn latex_to_proxy(tex: &str) -> String {
    let mut s = tex.to_string();
    for (from, to) in [
        (r"\mathrm", ""),
        (r"\textrm", ""),
        (r"\textsc", ""),
        (r"\rm", ""),
        (r"\tau", "τ"),
        (r"\rho", "ρ"),
        (r"\sigma", "σ"),
        (r"\odot", "⊙"),
        (r"\log", "log"),
        (r"\,", ""),
        (r"\;", ""),
        (r"\ ", " "),
        ("$", ""),
        ("{", ""),
        ("}", ""),
        ("_", ""),
        ("^", ""),
    ] {
        s = s.replace(from, to);
    }
    s
}

pub fn draw_legend<DB: DrawingBackend>(
    plot_area: &DrawingArea<DB, Shift>, // chart.plotting_area().strip_coord_spec()
    entries: &[LegendEntry<DB>],
    font: FontDesc<'static>,
    anchor: LegendAnchor,
    margin: i32,
    width_override: Option<i32>, // your manual tuning knob
) -> Result<(), Box<dyn Error>>
where
    DB::ErrorType: 'static,
{
    const PAD: i32 = 8;
    const MARKER_W: i32 = 24;
    const MARKER_TEXT_GAP: i32 = 8;
    const ROW_GAP: i32 = 6;

    let font_h = font.get_size() as i32;
    let row_h = font_h + ROW_GAP;
    // VAnchor::Top keeps inkscape's pdf_tex in the [t]-tabular construction
    // we've already calibrated — do NOT use VAnchor::Center here.
    let text_style = font
        .into_text_style(plot_area)
        .pos(Pos::new(HPos::Left, VPos::Top));

    let mut label_w: i32 = 0;
    for e in entries {
        let (w, _) = plot_area.estimate_text_size(&e.proxy, &text_style)?;
        label_w = label_w.max(w as i32);
    }

    let box_w = width_override.unwrap_or(PAD + MARKER_W + MARKER_TEXT_GAP + label_w + PAD);
    let box_h = PAD + entries.len() as i32 * row_h - ROW_GAP + PAD;

    let (pw, _) = plot_area.dim_in_pixel();
    let (x0, y0) = match anchor {
        LegendAnchor::UpperRight => (pw as i32 - margin - box_w, margin),
        LegendAnchor::UpperLeft => (margin, margin),
    };
    let (x1, y1) = (x0 + box_w, y0 + box_h);

    // Same look as configure_series_labels: translucent white bg, black border
    plot_area.draw(&Rectangle::new(
        [(x0, y0), (x1, y1)],
        WHITE.mix(0.8).filled(),
    ))?;
    plot_area.draw(&Rectangle::new([(x0, y0), (x1, y1)], BLACK.stroke_width(1)))?;

    let text_x = x0 + PAD + MARKER_W + MARKER_TEXT_GAP;
    for (i, e) in entries.iter().enumerate() {
        let row_top = y0 + PAD + i as i32 * row_h;
        let marker_y = row_top + font_h / 2;

        match &e.marker {
            LegendMarker::Line(c) => plot_area.draw(&PathElement::new(
                vec![(x0 + PAD, marker_y), (x0 + PAD + MARKER_W, marker_y)],
                c.stroke_width(2),
            ))?,
            LegendMarker::Dashed(c, on, off) => {
                let mut x = x0 + PAD;
                let x_end = x0 + PAD + MARKER_W;
                while x < x_end {
                    let xe = (x + on).min(x_end);
                    plot_area.draw(&PathElement::new(
                        vec![(x, marker_y), (xe, marker_y)],
                        c.stroke_width(2),
                    ))?;
                    x += on + off;
                }
            }
            LegendMarker::Custom(f) => f(plot_area, marker_y, x0 + PAD, x0 + PAD + MARKER_W)?,
        }

        plot_area.draw(&Text::new(
            e.tex.clone(),
            (text_x, row_top),
            text_style.clone(),
        ))?;
    }
    Ok(())
}
