#!/bin/bash

# Проверка наличия аргумента командной строки
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <input_file>"
    exit 1
fi

# Имя входного файла из аргумента командной строки
input_file="$1"

# Проверка существования входного файла
if [ ! -f "$input_file" ]; then
    echo "Error: File '$input_file' not found!"
    exit 1
fi

# Имя выходного файла на основе входного файла
output_file="${input_file%.*}.bed"

# Пропустить строки до [Probes]
skip_lines=true

# Открыть файл для записи
exec 3>"$output_file"

# Обработать строку файла
while IFS=$'\t' read -r Target_Region_Name Target_Region_ID Target_ID Species Build_ID Chromosome Start_Position End_Position Strand ULSO_Sequence ULSO_Genomic_Hits DLSO_Sequence DLSO_Genomic_Hits Probe_Strand Designer Design_Score Expected_Amplifed_Region_Size SNP_Masking; do
    # Проверить, достигли ли мы строки [Probes]
    if $skip_lines; then
        [[ "$Target_Region_Name" == "[Probes]" ]] && skip_lines=false
        continue
    fi
    
    # Записать строки в файл в формате BED
    echo -e "${Chromosome}\t${Start_Position}\t${End_Position}\t${Target_ID}\t${Design_Score}\t${Strand}" >&3
done < "$input_file"

# Закрыть файл
exec 3>&-

echo "Конвертация завершена. Файл сохранен как $output_file."

