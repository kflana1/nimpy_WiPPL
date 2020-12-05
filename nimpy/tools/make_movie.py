import argparse
import subprocess
from datetime import date
import os
"""-y overwrites file without asking.  I will be doing this."""


def main(folder, file_pattern, name, fps):
    if folder is None:
        folder = os.getcwd()
    else:
        folder = os.path.abspath(folder)

    if name is None:
        today = date.today()
        today_str = today.strftime("%Y-%m-%d")
        name = "{0:s}_movie.mp4".format(today_str)

    try:
        folder = os.path.normpath(folder)
        file_format = os.path.join(folder, file_pattern)
        call = ['ffmpeg', '-framerate', '{0}'.format(int(fps)), '-i', file_format, '-s', '1920x1080', '-c:v', 'libx264',
                '-pix_fmt', 'yuv420p', name, '-y']
        _ = subprocess.check_call(call)
    except subprocess.CalledProcessError as e:
        raise e

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Combine saved plots to a mp4 movie using ffmpeg.")
    parser.add_argument("--folder", "-f", type=str, help="Folder location of images to be combined to make movie.")
    parser.add_argument("--format", type=str, default="frame_%05d.png", help="File format string for ffmpeg")
    parser.add_argument("--name", '-n', type=str, help="Output movie filename (needs .mp4 at end)")
    parser.add_argument("--framerate", '-fps', type=int, default=30, help="Framerate for the movie")
    args = parser.parse_args()

    main(args.folder, args.format, args.name)
