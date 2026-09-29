from pathlib import Path

from PIL import Image


def main() -> None:
    source = Path("/scratch/td2248/projects/fixed_reference_credit_ge_20260929/ge_results/solve_v1/root_03/standard_diagnostics")
    destination = Path("/scratch/td2248/projects/fixed_reference_credit_ge_20260929/verification_v1/contact_sheets")
    destination.mkdir(parents=True, exist_ok=True)
    files = sorted(source.glob("*.png"))
    if len(files) != 17:
        raise RuntimeError(f"expected 17 standard PNGs, found {len(files)}")
    for sheet_index, start in enumerate(range(0, len(files), 4), start=1):
        images = [Image.open(path).convert("RGB") for path in files[start : start + 4]]
        width = max(image.width for image in images)
        height = max(image.height for image in images)
        canvas = Image.new("RGB", (2 * width, 2 * height), "white")
        for image_index, image in enumerate(images):
            canvas.paste(image, ((image_index % 2) * width, (image_index // 2) * height))
        canvas.save(destination / f"contact_sheet_{sheet_index:02d}.png")


if __name__ == "__main__":
    main()
