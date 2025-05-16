def write_file(file_name):
    alt_sfc = 0.0
    alt_top = 9270.0
    num_levels = 464.5
    with open(file_name, "w") as f:
        altitude = alt_sfc
        while altitude < alt_top:
            f.write(str(altitude) + '\n')
            altitude += (alt_top - alt_sfc)/(num_levels-1)
        f.write(str(alt_top))

def main():
    write_file('test.grd')

if __name__ == '__main__':
    main()