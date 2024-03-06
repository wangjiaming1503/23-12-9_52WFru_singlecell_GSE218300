install.packages("languageserver") 
.libPaths()           # 检查库路径
getOption("repos")    # 检查 CRAN 镜像地址
Sys.getenv("http_proxy")  # 检查 http 代理设置
Sys.getenv("https_proxy")  # 检查 https 代理设置
Sys.getenv("TZ")           # 检查时区设置
