input_file = 'ISAG.bed'
output_file = 'corrected_ISAG.bed'

with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    for line in infile:
        chrom, start, end = line.strip().split()
        start = int(start)  # Нулевое смещение не требуется, так как координаты для одиночных позиций
        end = int(end) + 1  # Конечная позиция должна быть увеличена на 1 для BED-формата
        outfile.write(f"{chrom}\t{start}\t{end}\n")

