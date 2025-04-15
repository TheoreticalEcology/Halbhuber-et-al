# copy Data, Results and Scripts
mastersThesis <- "/home/isabellehalbhuber/Toxicology/Scripts"
paper <- "/home/isabellehalbhuber/paper/EcoToxTM/Scripts"
file.copy(list.files(mastersThesis, full.names = TRUE), paper, recursive = TRUE)

