use std::error::Error;
use std::path::Path;
use std::process::Command;

const PNG_DPI: &str = "2000"; // rasterization density for PNGs

/// Full typeset-figure pipeline:
///   <stem>.svg ──inkscape --export-latex──▶ <stem>.pdf (slabs) + <stem>.pdf_tex
///   <stem>.pdf_tex ──patch──▶ corrected label anchoring (this is what the paper \inputs)
///   <stem>_standalone.tex ──pdflatex──▶ <stem>_standalone.pdf (self-contained, fonts embedded)
///
/// `font_px` = the base label size used in plotters (pt = px * 0.75 at 96 dpi).
pub fn svg_to_pdf(output_path: &str, font_px: f64) -> Result<(), Box<dyn Error>> {
    let svg_path = format!("{output_path}.svg");
    let pdf_path = format!("{output_path}.pdf");
    let pdf_tex_path = format!("{output_path}.pdf_tex");

    // 1. SVG -> slabs PDF + .pdf_tex
    let status = Command::new("inkscape")
        .args([
            "--export-type=pdf",
            "--export-latex",
            "--export-filename",
            &pdf_path,
            &svg_path,
        ])
        .status()?;
    if !status.success() {
        return Err(format!("inkscape failed on {svg_path}: {status}").into());
    }

    // 2. Fix label anchoring
    let n = patch_pdf_tex(Path::new(&pdf_tex_path))?;
    println!("patched {n} labels in {pdf_tex_path}");

    // 3. Standalone wrapper
    let path = Path::new(output_path);
    let dir = path
        .parent()
        .filter(|d| !d.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    let stem = path
        .file_name()
        .and_then(|s| s.to_str())
        .ok_or("svg_to_pdf: output_path has no usable file name")?;

    let font_pt = font_px * 0.75;
    let wrapper = format!(
        "\\documentclass{{standalone}}\n\
         \\usepackage{{graphicx, color}}\n\
         \\begin{{document}}%\n\
         \\sffamily\\fontsize{{{font_pt}}}{{{baseline:.1}}}\\selectfont%\n\
         \\input{{{stem}.pdf_tex}}%\n\
         \\end{{document}}\n",
        baseline = font_pt * 1.2,
    );
    let wrapper_name = format!("{stem}_standalone.tex");
    std::fs::write(dir.join(&wrapper_name), wrapper)?;

    // 4. Compile with cwd = figure dir so \input and the internal \includegraphics resolve
    let status = Command::new("pdflatex")
        .args(["-interaction=nonstopmode", "-halt-on-error", &wrapper_name])
        .current_dir(dir)
        .status()?;
    if !status.success() {
        return Err(format!(
            "pdflatex failed on {wrapper_name} — see {}",
            dir.join(format!("{stem}_standalone.log")).display()
        )
        .into());
    }

    // 5. Tidy aux files (keep the .tex — tiny, and useful for debugging)
    for ext in ["aux", "log"] {
        let _ = std::fs::remove_file(dir.join(format!("{stem}_standalone.{ext}")));
    }

    Ok(())
}

/// Rasterize the *typeset* PDF (not the SVG), so the PNG has the final
/// LaTeX-rendered labels. Runs the full pipeline first, so calling this
/// alone is sufficient.
pub fn svg_to_png(output_path: &str, font_px: f64) -> Result<(), Box<dyn Error>> {
    svg_to_pdf(output_path, font_px)?;

    let standalone_pdf = format!("{output_path}_standalone.pdf");
    let png_path = format!("{output_path}.png");

    let status = Command::new("inkscape")
        .args([
            &standalone_pdf,
            "--export-type=png",
            "--export-filename",
            &png_path,
            &format!("--export-dpi={PNG_DPI}"),
        ])
        .status()?;
    if status.success() {
        Ok(())
    } else {
        Err(format!("inkscape PNG export failed on {standalone_pdf}: {status}").into())
    }
}

/// Fix inkscape's --export-latex output in place: on label lines (the ones
/// containing a tabular), drop the exporter's `\lineheight{1.25}` and unwrap
/// `\smash{...}`, restoring the intended top-anchored geometry. Idempotent.
fn patch_pdf_tex(pdf_tex_path: &Path) -> Result<usize, Box<dyn Error>> {
    let content = std::fs::read_to_string(pdf_tex_path)?;
    let mut n = 0;

    let patched: String = content
        .split_inclusive('\n') // preserves line endings exactly
        .map(|line| {
            if line.contains(r"\begin{tabular}") {
                n += 1;
                strip_smash(&line.replace(r"\lineheight{1.25}", ""))
            } else {
                line.to_string()
            }
        })
        .collect();

    std::fs::write(pdf_tex_path, patched)?;
    Ok(n)
}

/// Remove every `\smash{...}` on a line, keeping the inner content.
/// Brace-matched, so nested `{...}` inside the content survives.
fn strip_smash(line: &str) -> String {
    const NEEDLE: &str = r"\smash{";
    let mut out = String::with_capacity(line.len());
    let mut rest = line;

    while let Some(start) = rest.find(NEEDLE) {
        out.push_str(&rest[..start]);
        let content_start = start + NEEDLE.len();

        let bytes = rest.as_bytes();
        let (mut depth, mut i) = (1usize, content_start);
        while i < bytes.len() && depth > 0 {
            match bytes[i] {
                b'{' => depth += 1,
                b'}' => depth -= 1,
                _ => {}
            }
            i += 1;
        }

        if depth == 0 {
            out.push_str(&rest[content_start..i - 1]); // inner content, sans closing brace
            rest = &rest[i..];
        } else {
            out.push_str(&rest[start..]); // unbalanced; leave untouched
            rest = "";
        }
    }
    out.push_str(rest);
    out
}
