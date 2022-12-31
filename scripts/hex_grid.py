import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection


def hex_grid(mat, cmap=None, file_name=None):
    shape = mat.shape
    i, j = np.meshgrid(range(shape[0]), range(shape[1]))
    i, j = i.flatten(), j.flatten()
    x, y = i + (j+1) % 2 / 2, j * np.sqrt(3) / 2
    centers = np.array([x, y]).T
    offsets = np.zeros((6, shape[0] * shape[1], 2))

    for i in range(6):
        offsets[i, :, 0] = np.sin(2 * np.pi * i / 6) / np.sqrt(3) * 1.01
        offsets[i, :, 1] = np.cos(2 * np.pi * i / 6) / np.sqrt(3) * 1.01
    verts = centers + offsets
    verts = np.swapaxes(verts, 0, 1)

    z = mat.flatten()

    fig, ax = plt.subplots()

    coll = PolyCollection(verts, array=z, cmap=cmap)
    ax.add_collection(coll)

    plt.axis('off')
    plt.autoscale(enable=True)
    fig.subplots_adjust(bottom=-0.06, top=1.06, left=-0.06, right=1.06)

    if file_name:
        plt.savefig(file_name, dpi=600)
        plt.close()
    else:
        plt.show()


def loadData(file_prefix, size=600):
    with open(f'{file_prefix}_vapor', 'rb') as file_stream:
        vapor_data = np.fromfile(
            file_stream, dtype=np.double).reshape(size, size)

    with open(f'{file_prefix}_solid', 'rb') as file_stream:
        solid_data = np.fromfile(
            file_stream, dtype=np.double).reshape(size, size)

    with open(f'{file_prefix}_liquid', 'rb') as file_stream:
        liquid_data = np.fromfile(
            file_stream, dtype=np.double).reshape(size, size)

    return vapor_data, solid_data, liquid_data


if __name__ == "__main__":
    inFolder = input('Input folder ("../dist/frames"): ')
    if len(inFolder) == 0:
        inFolder = '../dist/frames'

    outFolder = input('Output folder ("../dist/figures"): ')
    if len(outFolder) == 0:
        outFolder = '../dist/figures'

    frame = input('Show single frame ([N]/"end"/frame): ')
    if len(frame) == 0 or frame.lower() == 'n':
        start = int(input('Start position: '))
        frames = int(input('No. of frames: '))
        skip = int(input('Skip count: '))
        # Save all frames
        for i in range(start, frames, skip):
            print(f'Processing frame: {i}', end='\r')
            v, s, l = loadData(f'{inFolder}/snowflake_{i}')
            hex_grid(v + s + l, cmap='Blues',
                     file_name=f'{outFolder}/snowflake_{i:06}.png')
    elif frame.lower() == 'end':
        # Show final image
        v, s, l = loadData(f'{inFolder}/snowflake')
        hex_grid(v + s + l, cmap='Blues')
    else:
        # Show specific frame
        frame = int(frame)
        v, s, l = loadData(f'{inFolder}/snowflake_{frame}')
        hex_grid(v + s + l, cmap='Blues')
