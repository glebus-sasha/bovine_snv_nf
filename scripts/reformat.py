import argparse
import re

# Словарь для преобразования хромосом
chromosome_mapping = {
    '1': 'AC_000158.1',
    '2': 'AC_000159.1',
    '3': 'AC_000160.1',
    '4': 'AC_000161.1',
    '5': 'AC_000162.1',
    '6': 'AC_000163.1',
    '7': 'AC_000164.1',
    '8': 'AC_000165.1',
    '9': 'AC_000166.1',
    '10': 'AC_000167.1',
    '11': 'AC_000168.1',
    '12': 'AC_000169.1',
    '13': 'AC_000170.1',
    '14': 'AC_000171.1',
    '15': 'AC_000172.1',
    '16': 'AC_000173.1',
    '17': 'AC_000174.1',
    '18': 'AC_000175.1',
    '19': 'AC_000176.1',
    '20': 'AC_000177.1',
    '21': 'AC_000178.1',
    '22': 'AC_000179.1',
    '23': 'AC_000180.1',
    '24': 'AC_000181.1',
    '25': 'AC_000182.1',
    '26': 'AC_000183.1',
    '27': 'AC_000184.1',
    '28': 'AC_000185.1',
    '29': 'AC_000186.1',
    'X': 'AC_000187.1'
}

# Функция для преобразования хромосом в строке BED файла
def convert_chromosome_format(bed_line, mapping):
    fields = bed_line.strip().split('\t')
    chrom_id = fields[0]
    if chrom_id in mapping:
        fields[0] = mapping[chrom_id]
    return '\t'.join(fields)

# Основная функция для обработки BED файла
def process_bed_file(input_bed, output_bed, mapping):
    with open(input_bed, 'r') as infile, open(output_bed, 'w') as outfile:
        for line in infile:
            converted_line = convert_chromosome_format(line, mapping)
            outfile.write(converted_line + '\n')

# Настройка аргументов командной строки
def main():
    parser = argparse.ArgumentParser(description='Convert chromosome names in a BED file.')
    parser.add_argument('input_bed', help='Input BED file')
    parser.add_argument('output_bed', help='Output BED file')
    args = parser.parse_args()
    
    process_bed_file(args.input_bed, args.output_bed, chromosome_mapping)

if __name__ == '__main__':
    main()

