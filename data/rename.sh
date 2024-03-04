for dir in ~/data_linli/bioinfo/2023-12-9_52WFru_singlecell_GSE218300/data/*/; do
    cd "$dir"
    # 对于WD前缀的文件
    for file in GSM*_WD_*; do
        mv "$file" "${file#GSM*_WD_}"
    done
    # 对于Chow前缀的文件
    for file in GSM*_Chow_*; do
        mv "$file" "${file#GSM*_Chow_*_}"
    done
    cd - # 返回之前的目录
done