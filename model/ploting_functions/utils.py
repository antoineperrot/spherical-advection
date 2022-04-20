from tqdm import tqdm 
import os
from IPython.display import clear_output

def save_plots(save_directory,img_name, plotting_function,times, titles, plotting_func_args):
     #cleaning directory :
    for file in os.listdir(save_directory):
        os.remove(save_directory + file)

    for index in tqdm(range(len(times)),position=0,leave=True) :
        fig = plotting_function(times[index], *plotting_func_args)
        fig.savefig(save_directory + img_name + str(index_img(index)) )
        fig.clear()

    clear_output(wait=False)
    
import glob
from PIL import Image


def make_gif(source_imgs_path, output_path, title):
    fp_in =source_imgs_path +'*.png'
    fp_out = output_path + title +'.gif'
    img, *imgs = [Image.open(f) for f in sorted(glob.glob(fp_in))]
    img.save(fp=fp_out, format='GIF', append_images=imgs,
             save_all=True, duration=120, loop=0)


def index_img(index):
    if index < 10 :
        return '00'+str(index) 
    elif 10<= index < 100 :
        return '0'+str(index) 
    else:
        return index
        
def make_video(imgs_dir, output_dir, output_filename, fps= 25):
    import cv2
    import glob

    img_array = []
    for filename in glob.glob(imgs_dir + '*.png'):
        img = cv2.imread(filename)
        height, width, layers = img.shape
        size = (width,height)
        img_array.append(img)

    out = cv2.VideoWriter( output_dir + output_filename + '.mp4', 		
    cv2.VideoWriter_fourcc(*'DIVX'), fps, size)

    for i in range(len(img_array)):
        out.write(img_array[i])
    out.release()

import numpy as np
def get_stop_times(times, n_stops):
    stops = np.linspace(0,len(times)-1,n_stops).astype(int)
    times_stops = [times[stop] for stop in stops]
    return times_stops 
        
