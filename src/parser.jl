function read_files(folder_name)

    println(folder_name)

    files_ = readdir(folder_name)
    # println(files_)

    files = [f for f in files_ if f[end-3:end] == ".cor"]

    println(files)

end

read_files("/home/mathieu/Documents/TSSP/tssp")