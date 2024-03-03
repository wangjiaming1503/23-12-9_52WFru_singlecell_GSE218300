# 创建示例数据框
df <- data.frame(original_codes = c("01", "11", "01", "11", "01"))

# 定义转换函数
convert_codes_to_labels <- function(df, column_name, new_column_name) {
  # 使用ifelse函数根据条件增加新列
  df[[new_column_name]] <- ifelse(df[[column_name]] == "01", "A", 
                                  ifelse(df[[column_name]] == "11", "B", NA))
  return(df)
}

# 应用转换函数
df <- convert_codes_to_labels(df, "original_codes", "new_labels")

# 查看转换后的数据框
print(df)
