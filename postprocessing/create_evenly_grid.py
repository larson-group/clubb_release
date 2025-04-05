def write_file(file_name):
    alt_sfc = 20.0
    alt_top = 9000.0
    num_levels = 33
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