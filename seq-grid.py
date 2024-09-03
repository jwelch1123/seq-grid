from PIL import Image, ImageDraw

def find_square_factors(n):
    '''
        Returns the factors of n that are closest to being a square.
        Straight up written by Copilot
    '''
    hold_diff: int = n
    hold_lower: int = None
    for i in range(1, n):
        if n % i == 0 and abs(int((n/i)-i)) <= hold_diff:
            hold_lower = i
            hold_diff = abs(int((n/i)-i))
    if abs(int(n/hold_lower)) == 1 or hold_lower == 1:
        print("This might be a prime")
    return abs(int(n/hold_lower)), hold_lower

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

def assemble_key(dims: tuple, seq: list, **kwargs):
    '''
        package sequence into a dictionary with the keys being the coordinates of the grid
    '''
    grid_dict = {}
    for i, val in enumerate(seq):
        x = i % dims[0]
        y = i // dims[0]
        grid_dict[(x, y)] = {"seq_val": val, **kwargs}
    return grid_dict

def add_spacing(grid_dict: dict, grid_dims: tuple, img_dims: tuple):
    '''
        Add spacing to the grid
    '''

    grid_x, grid_y = grid_dims
    img_x, img_y = img_dims

    #print("grid Dims: ", grid_dims)
    #print("Img Dims: ", img_dims)
    x_spacing = int(img_x / (2 + grid_x))
    y_spacing = int(img_y / (2 + grid_y))

    #print("Spacing: ", x_spacing, y_spacing)

    for key in grid_dict:
        x1 = (key[0]+1) * x_spacing
        y1 = (key[1]+1) * y_spacing

        grid_dict[key]["coords"] = (x1, y1)
    
    return grid_dict

def bottom_up(dims: tuple, seq: list, shape_size: int, shape_pad: int, ext_margin: int):
    '''
        Fully assembly of the shape grid in 1 pass without lookups. 
    '''

    grid_dict = {}
    for i, val in enumerate(seq):
        x = i % dims[0]
        y = i // dims[0]
        
        x_coord = ext_margin + (x * (shape_size + shape_pad))
        y_coord = ext_margin + (y * (shape_size + shape_pad))

        grid_dict[(x, y)] = {"seq_val": val, 
                             "shape_size": shape_size, 
                             "coords":(x_coord, y_coord)
                             }

    return grid_dict

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

def draw_grid(grid_dict: dict, dims: tuple, color_dict: dict, bkgd_color: str = "black"):
    '''
        Draws the grid
    '''

    img = Image.new('RGB', dims, color = bkgd_color) 
    draw = ImageDraw.Draw(img)

    for key in grid_dict:
        x1 = grid_dict[key]['coords'][0]
        y1 = grid_dict[key]['coords'][1]
        x2 = x1 + 10
        y2 = y1 + 10
        color = color_dict[grid_dict[key]["seq_val"]]
        #print(f"Coords calc: grid_dict[key]: {grid_dict[key]}, coords: {grid_dict[key]['coords']}")
        #print("Drawing: ", x1, y1, x2, y2)

        #draw.rectangle([x1, y1, x2, y2], fill = color)
        draw.circle(((x1+x2)/2, (y1+y2)/2),radius=(x2-x1)/2, fill = color)
    
    img.show()

    return img

def main():

    seq = "ATGCGTACGATCGTAG"
    seq = "ATGTTCTCTCCAATTTTGTCCTTGGAAATTATTTTAGCTTTGGCTACTTTGCAATCTGTC"

    #seq = get_genbank_sequence("NC_000913.3")
    #seq = read_txt("Test Sequences\collagen.txt")
    #seq = read_txt("Test Sequences\CoVid-19.txt")
    seq = read_txt("Test Sequences\Truncated_100k_refChromo1.txt")
    #seq = read_txt("Test Sequences/NC_000001.11[77458449..201936659].fa")

    print("*"*50 + "\n" + f"Sequence length: {len(seq)}" + "\n")
    
    # image ratio
    x_ratio = 1
    y_ratio = 1

    shape_size = 10 
    shape_pad = 2
    ext_margin = 10


    x, y = find_ratio_factors(len(seq), x_ratio, y_ratio)
    #x, y = find_square_factors(len(seq))
    
    print("*"*50 + "\n" + f"X, Y: {x},{y}")

    # image dimensions, something is wrong here, Padding added to top but not bottom?
    img_x =  (x * shape_size) + ((x-1) * shape_pad) + (2 * ext_margin)
    img_y =  (y * shape_size) + ((y-1) * shape_pad) + (2 * ext_margin)

    print("*"*50 + "\n" + f"Img Dims: {img_x},{img_y}")

    #prot_color_dict = {}
    color_dict = {"A": "Green", "T": "Red", "C": "Orange", "G": "Blue"}
   
    #import time
    #print("Time:", time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

    spaced_grid_dict = bottom_up((x,y), seq, shape_size, shape_pad, ext_margin)

    #print("*"*50 + "\n" + f"Spaced key: {spaced_grid_dict[(0,0)]}")
    #print("Time:", time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    
    img = draw_grid(spaced_grid_dict, (img_x, img_y), color_dict)
    
    #print("Time:", time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    
    #img.save("./in progress/____.png")

    #img.show()
    pass


if __name__ == "__main__":
    main()