import os
from urllib.request import urlopen

from PIL import Image
from PIL import ImageFont
from PIL import ImageDraw
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

from docu.visualizations import plot_rpolygon_max_rectangles
from rportion import rempty, rclosed


def main():

    for filename, text in [("logo_single.png", "R"), ("logo_double.png", "RP")]:

        a = urlopen(
            "https://github.com/google/fonts/blob/90d7886db9000c893b9559828bf028aaed5f9c10/ofl/bungee/Bungee-Regular.ttf?raw=true")
        font = ImageFont.truetype(a, 30)
        _, _, w, h = font.getbbox(text)
        h *= 2
        image = Image.new('L', (w, h), 1)
        draw = ImageDraw.Draw(image)
        draw.text((0, 0), text, font=font)
        char_image = np.asarray(image)
        char_image = np.where(char_image, 0, 1)
        char_image = char_image[(char_image != 0).any(axis=1)]

        downscaled = char_image[::1, ::1]

        plt.imshow(char_image)
        plt.show()
        plt.imshow(downscaled)
        plt.show()

        print("Create Polygon")
        poly = rempty()
        for i in tqdm(range(downscaled.shape[0])):
            for j in range(downscaled.shape[1]):
                if downscaled[i, j]:
                    poly |= rclosed(j, j+1, i, i+1)

        print("Plot Polygon")
        fig, ax = plt.subplots(1, 1, figsize=(15*len(text), 15))
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.get_xaxis().set_ticks([])
        ax.get_yaxis().set_ticks([])
        plot_rpolygon_max_rectangles(ax, poly, box=None, alpha_used=1.0, plot_boundary=True, verbose=True)
        ax.invert_yaxis()
        fig.savefig(filename, bbox_inches='tight', pad_inches=0)


if __name__ == "__main__":
    main()
