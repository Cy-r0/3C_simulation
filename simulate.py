#!/usr/bin/env python3
"""Chromosome Conformation Capture (3C) simulation"""


__author__ = 'Ciro Cursio'


import math

from PIL import Image, ImageDraw
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sn


class Simulator():

    def __init__(self, length, angle_variance):
        self.filament = self.generate_filament(length, angle_variance)
        self.overlap_map = self.calc_overlaps()

    def generate_filament(self, l=10000, angle_sigma=0.15):
        """Generate a sequence of 2D points with these constraints:
        1. Each pair of consecutive points is connected (i.e. there's a small
            constant distance between them);
        2. The direction of the vector between consecutive points is mostly random,
            but not completely: to encourage smooth curves, the direction between
            each pair is highly correlated with the direction between neighbouring
            pairs.
        Note, position is quantised on a grid and angle is quantised to Von Neumann
            neighbourhood.
        Args:
            l: filament length.
            angle_sigma: variance of change in angle, expressed as fraction of 2pi
                (higher sigma = angle changes more sharply).
        """
        print('Generating filament...')
        filament = np.zeros((l, 2), dtype='int')
        angle = 2 * math.pi * np.random.random()
        for i in range(1, l):
            angle_q = (angle + math.pi/4) // (math.pi/2) * math.pi/2
            v = [int(math.cos(angle_q)), int(math.sin(angle_q))]
            filament[i] = filament[i-1] + v
            angle += (np.random.random() - 0.5) * 2 * math.pi * angle_sigma
        print('Done.')
        return filament

    def calc_overlaps(self):
        """Calculate which points of the filament overlap."""
        print('Calculating overlaps...')
        overlap_map = []
        for i in range(len(self.filament)):
            point = self.filament[i]
            all_i = np.squeeze(np.argwhere(np.all(self.filament==point, axis=-1)), axis=-1)
            # Remove idx of current point
            all_i = np.delete(all_i, np.where(all_i==i))
            overlap_map.append(all_i)
        print('Done.')
        return overlap_map

    def _scale_filament(self, filament, img_size, padding):
        """
        Args:
            filament: chromatin array.
            img_size: size of image.
            padding: pixels around filament bounds.
        """
        filament = self.filament.copy().astype('float')
        min_x = np.min(filament[:, 0])
        max_x = np.max(filament[:, 0])
        min_y = np.min(filament[:, 1])
        max_y = np.max(filament[:, 1])

        filament[:, 0] = (filament[:, 0] - min_x) / (max_x - min_x)
        filament[:, 1] = (filament[:, 1] - min_y) / (max_y - min_y)
        if (max_x - min_x) > (max_y - min_y):
            filament[:, 1] *= (max_y - min_y) / (max_x - min_x)
        else:
            filament[:, 0] *= (max_x - min_x) / (max_y - min_y)

        return filament * (img_size - 2 * padding) + padding

    def visualise_filament(self, img_size=1000, padding=50, highlight_intersections=False):
        filament = self._scale_filament(self.filament, img_size, padding)
        img = Image.new('RGB', (img_size, img_size), '#303030')
        overlay = Image.new('RGBA', img.size, '#00000000')
        draw = ImageDraw.Draw(img)
        draw_overlay = ImageDraw.Draw(overlay)
        color_start = (255, 128, 0)
        color_end = (255, 0, 0)
        for i in range(len(filament)):
            x = filament[i, 0]
            y = filament[i, 1]
            if len(self.overlap_map[i]) and highlight_intersections:
                draw_overlay.rectangle((x-2, y-2, x+2, y+2), fill='cyan')
            else:
                a = i/len(filament)
                color = (int((1 - a) * color_start[0] + a * color_end[0]),
                        int((1 - a) * color_start[1] + a * color_end[1]),
                        int((1 - a) * color_start[2] + a * color_end[2]))
                draw.rectangle((x-1, y-1, x+1, y+1), fill=color)
        img.paste(overlay, overlay)
        img.show()
        img.save('filament.png')

    def simulate_3C(self, segment_l):
        n = len(self.filament) // segment_l + 1
        matrix = np.zeros((n, n))
        for i in range(len(self.filament)):
            segment_i = i // segment_l
            if len(self.overlap_map[i]):
                for j in self.overlap_map[i]:
                    # calculate index of segment the overlap falls in
                    segment_j = j // segment_l
                    matrix[segment_i, segment_j] = 1
        np.fill_diagonal(matrix, 2)
        self.contact_matrix = matrix
        matrix_plot = sn.heatmap(self.contact_matrix, cmap=sn.cm.rocket_r)
        plt.show()
        fig = matrix_plot.get_figure()
        fig.savefig('contact_matrix.png')


if __name__ == '__main__':
    # Good 10k long filaments:
    # seed 2, variance 0.2 - two nests separated by single strand
    # seed 9, variance 0.1 - a few short range interactions, but not much else
    # seed 4, variance 0.3 - one medium range interaction
    # seed 10, variance 0.3 - long range interaction
    np.random.seed(10)
    sim = Simulator(length=10000, angle_variance=0.3)
    sim.visualise_filament(highlight_intersections=True)
    sim.simulate_3C(segment_l=100)
