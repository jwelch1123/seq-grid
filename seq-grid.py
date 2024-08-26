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

def get_genbank_sequence(accession: str):
    pass

def draw_grid(grid_dict: dict, dims: tuple, bk_color: str = "black"):
    '''
        Draws the grid
    '''
    x, y = dims

    img = Image.new('RGB', (x*100, y*100), color = bk_color) 
    draw = ImageDraw.Draw(img)

    for key in grid_dict:
        x1 = key[0]
        y1 = key[1]
        x2 = x1
        y2 = y1
        draw.rectangle([x1, y1, x2, y2], fill = grid_dict[key]["color"])
    
    img.show()

    return img

def main():
    # image dimensions
    img_x = 1000
    img_y = 1000
    
    # image ratio
    x_ratio = 1
    y_ratio = 1



    #seq = get_genbank_sequence("NC_000913.3")
    #DNA_color_dict = {"A": "#FF0000", "T": "#FFFF00", "C": "#0000FF", "G": "#008000"}
    #prot_color_dict = {}
    color_dict = {"A": "Green", "T": "Red", "C": "Orange", "G": "Blue"}
    seq = "ATGCGTACGATCGTAG"
    seq = "ATGTTCTCTCCAATTTTGTCCTTGGAAATTATTTTAGCTTTGGCTACTTTGCAATCTGTC"

    n = len(seq)
    
    x, y = find_ratio_factors(n, x_ratio, y_ratio)
    print(f"X: {x}, Y: {y}")
    
    grid_dict = assemble_key((x,y), seq)
    print(f"Grid Dict: {grid_dict}\n\n")

    colored_grid_dict = apply_colors(grid_dict, color_dict)
    print(f"Colored Grid Dict: {colored_grid_dict}")
    


    img = draw_grid(colored_grid_dict, (img_x, img_y))

    pass


if __name__ == "__main__":
    main()