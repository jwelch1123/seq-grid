from PIL import Image, ImageDraw
import time

def find_ratio_factors(n: int, x:int, y:int):
    '''
        Returns the scale factor for x & y of the smallest grid large enough to contain n items 
    '''
    scale_factor = 1
    while (x * y * (scale_factor ** 2)) < n:
        scale_factor += 1
        if scale_factor > n:
            break

    return x * scale_factor, y * scale_factor 

def get_genbank_sequence(accession: str):
    pass

def read_txt(file_path: str):
    if file_path.endswith(".fa"):
        with open(file_path, "r") as file:
            lines = file.readlines()
            seq = ''.join(lines[1:]).strip()
            seq = ''.join(filter(str.isalpha, seq))

    else:
        
        with open(file_path, "r") as file:
            seq = file.read().strip()
            seq = ''.join(filter(str.isalpha, seq))
    return seq

def draw_grid(dims: tuple, seq: list, shape_size: int, shape_pad: int, ext_margin: int, img_dims:tuple, color_dict: dict, bkgd_color: str = "black"):
    '''
        Draws the grid
    '''
    img = Image.new('RGB', img_dims, color = bkgd_color) 
    draw = ImageDraw.Draw(img)

    for i, val in enumerate(seq):

        if i % 1000000 == 0:
            print(f"Hit: {i/1000000}M")
            print("at Time:", time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

        x = i % dims[0]
        y = i // dims[0]
        
        x_coord = ext_margin + (x * (shape_size + shape_pad))
        y_coord = ext_margin + (y * (shape_size + shape_pad))

        x2 = x_coord + shape_size
        y2 = y_coord + shape_size

        color = color_dict[val]
       
        #draw.rectangle([x_coord, y_coord, x2, y2], fill = color)
        draw.circle(((x_coord+x2)/2, (y_coord+y2)/2),radius=(x2-x_coord)/2, fill = color)
    
    img.show()

    return img

def main():

    #seq = "ATGCGTACGATCGTAG"
    #seq = "ATGTTCTCTCCAATTTTGTCCTTGGAAATTATTTTAGCTTTGGCTACTTTGCAATCTGTC"
    #seq = get_genbank_sequence("NC_000913.3")
    #seq = read_txt("Test Sequences\collagen.txt")
    #seq = read_txt("Test Sequences\CoVid-19.txt")
    #seq = read_txt("Test Sequences\Truncated_100k_refChromo1.txt")
    #seq = read_txt("Test Sequences\Truncated_10M_refChromo1.txt")
    seq = read_txt("Test Sequences/NC_000001.11[77458449..201936659].fa")

    print("*"*50 + "\n" + f"Sequence length: {len(seq)}" + "\n")

    # image ratio
    x_ratio = 1
    y_ratio = 1

    shape_size = 10 
    shape_pad = 2
    ext_margin = 10

    x, y = find_ratio_factors(len(seq), x_ratio, y_ratio)
    #x, y = find_square_factors(len(seq))

    img_x =  (x * shape_size) + ((x-1) * shape_pad) + (2 * ext_margin)
    img_y =  (y * shape_size) + ((y-1) * shape_pad) + (2 * ext_margin)


    #prot_color_dict = {}
    color_dict = {"A": "Green", "T": "Red", "C": "Orange", "G": "Blue"}
   

    print("Start Time:", time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    
    img = draw_grid((x,y), seq, shape_size, shape_pad, ext_margin, (img_x, img_y), color_dict) # forgot to implement shape size

    print("End Time:", time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

    #print("Time:", time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    
    #img.save("./in progress/human_ref_chromosome1.png")

    #img.show()
    pass


if __name__ == "__main__":
    main()