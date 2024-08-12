import pandas as pd

# Чтение файла corrected_output.bed
bed_df = pd.read_csv("corrected_output.bed", sep='\t', header=None, names=["chromosome", "start", "end", "name"])

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

# Удаление 'chr' из названий хромосом
bed_df['chromosome'] = bed_df['chromosome'].str.replace("chr", "")

# Замена значений хромосом по словарю
bed_df['chromosome'] = bed_df['chromosome'].map(chromosome_mapping)

# Сохранение результата в новый BED файл
bed_df.to_csv("final_remapped_output.bed", sep='\t', index=False, header=False)

