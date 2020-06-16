module peptideutils

    export translate1to3, translate3to1

    peptideutils_letters1 = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
    peptideutils_letters3 = ["ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HSE", "ILE", "LYS", "LEU", "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR"]

    function translate1to3(string)
        code = split(string, "")
        new_string = ""
        for letter in code
            index = findall(x->x==letter, peptideutils_letters1)[1]
            new_string = new_string * peptideutils_letters3[index] * "-"
        end
        #println(new_string)
        str_length = length(new_string)
        new_string = new_string[1:str_length-1]
        return new_string
    end

    function translate3to1(string)
        code = split(string, "-")
        new_string = ""
        for AA in code
            if AA == "HIS"
                AA = "HSE"
            end
            index = findall(x->x==AA, peptideutils_letters3)[1]
            new_string = new_string * peptideutils_letters1[index]
        end
        return new_string
    end

end
