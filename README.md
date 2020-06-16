# Judred

Judred is a peptide specific molecular descriptor generator, it is loosely based on the [Mordred](https://github.com/mordred-descriptor/mordred) python package. The remit of this programme is to generate descriptors of peptides from only single-letter codes in order to provide a significantly faster method of generating descriptors for very large search spaces. It uses the HDF5 file format to store the data as this format is supported by many programming languages and platforms.

## Example usage

Generate the full dataset of zwitterionic dipeptides:
```julia
julia judred.jl 2 
```

Generate the full dataset of zwitterionic tripeptides:
```julia
julia judred.jl 3
```

Generate the full dataset of zwitterionic tetrapeptides:
```julia
julia judred.jl 4 
```

etc.

## Example output
![HDFView dipeptides](https://raw.githubusercontent.com/avanteijlingen/Judred/master/HDFView.PNG)