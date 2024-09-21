import os
from dotenv import load_dotenv, find_dotenv
from PIL import Image, ImageDraw
from Bio import Entrez, SeqIO
import math

load_dotenv(find_dotenv())

def find_ratio_factors(l: int, x:int, y:int, tight: bool):
    '''
        Returns the scale factor for x & y of the smallest grid large enough to contain n items 
    '''
    scale_factor = 1
    while (x * y * (scale_factor ** 2)) < l:
        scale_factor += 1
        if scale_factor > l:
            break

    x = x * scale_factor
    y = y * scale_factor

    if tight:
        y = math.ceil(l / x)
        x = l if l < x else x
    
    return x, y 

def get_genbank_sequence(accession: str):
    Entrez.email = os.getenv("entrez_email")
    handle = Entrez.efetch(
                        db="nucleotide", 
                        id=accession, 
                        rettype="gb", 
                        retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()
    return record.seq

def read_txt(file_path: str):
    if file_path.endswith(".fa"):
        with open(file_path, "r") as file:
            lines = file.readlines()
            seq = ''.join(lines[1:]).strip()
            seq = ''.join(filter(str.isalpha, seq))

    elif file_path.endswith(".txt"):
        
        with open(file_path, "r") as file:
            seq = file.read().strip()
            seq = ''.join(filter(str.isalpha, seq))
    return seq

def validate_seq(seq: str, valid_chars: list):    
    return all([char.lower() in valid_chars for char in seq])

def length_control(seq: str, max_len: int):
        return len(seq) > max_len

def draw_grid(dims: tuple, seq: list, shape_size: int, shape_pad: int, ext_margin: int, color_dict: dict, bkgd_color: str = "black"):
    '''
        Draws the grid
    '''

    img_width  = (ext_margin * 2) + (dims[0] * shape_size) + ((dims[0] - 1) * shape_pad) #- shape_pad 
    img_height = (ext_margin * 2) + (dims[1] * shape_size) + ((dims[1] - 1) * shape_pad) #- shape_pad
    
    print("Image Dimensions: ", img_width, img_height)

    img = Image.new('RGB', (img_width, img_height), color = bkgd_color) 
    draw = ImageDraw.Draw(img)

    for i, val in enumerate(seq):

        x = i % dims[0] 
        y = i // dims[0]
        
        x_coord = ext_margin + (x * (shape_size + shape_pad))
        y_coord = ext_margin + (y * (shape_size + shape_pad))

        x2 = x_coord + shape_size-1
        y2 = y_coord + shape_size-1

        color = color_dict[val]

        print("For val: ", val)
        print("x,y: ", x, y)
        print("coords: ", x_coord, y_coord)
        print("x2,y2: ", x2, y2)       
        print("margin, pad: ", ext_margin, shape_pad)
        print("\n")
        #draw.rectangle([x_coord, y_coord, x2, y2], fill = color)
        draw.circle(((x_coord+x2)/2, (y_coord+y2)/2),radius=(x2-x_coord)/2, fill = color)

    return img

def seq_generator(seq):
    for item in seq:
        yield item

def main():

    #seq = "ATGCGTACGATCGTAG"
    #seq = "ATGTTCTCTCCAATTTTGTCCTTGGAAATTATTTTAGCTTTGGCTACTTTGCAATCTGTC"
    seq = 'ACGTACGTACGTACGT'
    #seq = get_genbank_sequence("NM_001301717")
    #seq = read_txt("Test Sequences\collagen.txt")
    #seq = read_txt("Test Sequences\CoVid-19.txt")
    #seq = read_txt("Test Sequences\Truncated_100k_refChromo1.txt")
    #seq = read_txt("Test Sequences\Truncated_10M_refChromo1.txt")
    #seq = read_txt("Test Sequences/dystrophin.fa")

    # Sequence Settings
    max_seq_len = 2500000 # 2.5M
    color_dict = {"A": "Green", "T": "Red", "C": "Orange", "G": "Blue"}
    #prot_color_dict = {} 
    # https://biosci.mcdb.ucsb.edu/biochemistry/tutorials/pdbtutorial/console/predefinedcolorsdiagram.html

    # image settings
    x_ratio = 1
    y_ratio = 1
    ext_margin = 1
    tight_img = False

    # shape settings
    shape_size = 2#10 
    shape_pad = 1
   


    if not validate_seq(seq, [key.lower() for key in color_dict.keys()]):
        ValueError("Sequence contains invalid characters, please check the sequence and try again")
    if not length_control(seq, max_seq_len):
        ValueError(f"Provided sequence is longer than allowed, please keep sequences to <{max_seq_len} characters")
    

    x, y = find_ratio_factors(len(seq), x_ratio, y_ratio, tight_img)

    seq_gen = seq_generator(seq)
    img = draw_grid((x,y), seq_gen, shape_size, shape_pad, ext_margin, color_dict)
    
    #img.save("./in progress/__.png")

    img.show()
    #print(list(img.getdata()))
    

if __name__ == "__main__":
    main()