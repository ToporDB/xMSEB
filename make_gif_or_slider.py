import sys
from pathlib import Path
from PIL import Image

def make_gif(source_dir, output_dir, gif_name="animation.gif", duration=300):
    frames = []
    images = sorted(Path(source_dir).glob("*.png"))

    for img_path in images:
        frames.append(Image.open(img_path))

    if frames:
        output_path = Path(output_dir) / gif_name
        frames[0].save(output_path, format='GIF', append_images=frames[1:], save_all=True, duration=duration, loop=0)
        print(f"GIF saved to {output_path}")
    else:
        print("No images found to create GIF!")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: make_gif_or_slider.py <screenshot_dir> <output_dir>")
        sys.exit(1)

    screenshot_dir = sys.argv[1]
    output_dir = sys.argv[2]

    make_gif(screenshot_dir, output_dir)
