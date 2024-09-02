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

def assemble_key(dims: tuple, seq: list, seq_title: str = "seq_val"):
    '''
        package sequence into a dictionary with the keys being the coordinates of the grid
    '''
    dims = [(x, y) for y in range(dims[1]) for x in range(dims[0]) ]
    grid_dict = {dims[loc]: {seq_title: value} for loc, value in enumerate(seq)}
    return grid_dict

def apply_colors(grid_dict, color_dict: dict):
    '''
        Apply colors to the dictionary of the grid
    '''
    for key in grid_dict:
        grid_dict[key]["color"] = color_dict[grid_dict[key]["seq_val"]]
    return grid_dict

def add_spacing(grid_dict: dict, grid_dims: tuple, img_dims: tuple):
    '''
        Add spacing to the grid
    '''

    grid_x, grid_y = grid_dims
    img_x, img_y = img_dims

    print("grid Dims: ", grid_dims)
    print("Img Dims: ", img_dims)
    x_spacing = int(img_x / (2 + grid_x))
    y_spacing = int(img_y / (2 + grid_y))

    print("Spacing: ", x_spacing, y_spacing)

    for key in grid_dict:
        x1 = (key[0]+1) * x_spacing
        y1 = (key[1]+1) * y_spacing

        grid_dict[key]["coords"] = (x1, y1)
    
    return grid_dict

def get_genbank_sequence(accession: str):
    pass

def draw_grid(grid_dict: dict, dims: tuple, bkgd_color: str = "black"):
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
        print(f"Coords calc: grid_dict[key]: {grid_dict[key]}, coords: {grid_dict[key]['coords']}")
        #print("Drawing: ", x1, y1, x2, y2)
        #draw.rectangle([x1, y1, x2, y2], fill = grid_dict[key]["color"])
        draw.circle(((x1+x2)/2, (y1+y2)/2),radius=(x2-x1)/2, fill = grid_dict[key]["color"])
    
    img.show()

    return img

def main():
    # image dimensions
    # can derive these from the idea factors
    img_x = 100
    img_y = 100  
    
    # image ratio
    x_ratio = 1
    y_ratio = 1



    #seq = get_genbank_sequence("NC_000913.3")
    #DNA_color_dict = {"A": "#FF0000", "T": "#FFFF00", "C": "#0000FF", "G": "#008000"}
    #prot_color_dict = {}
    color_dict = {"A": "Green", "T": "Red", "C": "Orange", "G": "Blue"}
    seq = "ATGCGTACGATCGTAG"
    seq = "ATGTTCTCTCCAATTTTGTCCTTGGAAATTATTTTAGCTTTGGCTACTTTGCAATCTGTC"

    #seq = open("Test Sequences\collagen.txt", "r").read().strip()
    #seq = ''.join(filter(str.isalpha, seq))

    n = len(seq)
    
    x, y = find_ratio_factors(n, x_ratio, y_ratio)
    #x, y = find_square_factors(n)
    #print(f"X: {x}, Y: {y}")
    
    grid_dict = assemble_key((x,y), seq)
    #print(f"Grid Dict: {grid_dict}\n\n")

    colored_grid_dict = apply_colors(grid_dict, color_dict)
    #print(f"Colored Grid Dict: {colored_grid_dict}")
    
    spaced_grid_dict = add_spacing(colored_grid_dict, (x,y), (img_x, img_y))
    print("Spaced Grid Dict: ", spaced_grid_dict[(0,0)], spaced_grid_dict[(2,2)])

    draw_grid(spaced_grid_dict, (img_x, img_y))

    #img.show()
    pass


if __name__ == "__main__":
    main()