# 获取所有已安装的包的信息
installed_packages_local <- installed.packages(lib.loc = system_library_paths[3]) 
installed_packages_local_site <- installed.packages(lib.loc = system_library_paths[2]) 
installed_packages_local_user <- installed.packages(lib.loc = system_library_paths[1]) 

# 获取系统库的路径
system_library_paths <- .libPaths()

# 列出系统库中的包
system_packages <- installed_packages[installed_packages["/usr/local/lib/R/site-library"] %in% system_library_paths, ]

# 打印系统库中的包的名称
print(system_packages[, "Package"])
