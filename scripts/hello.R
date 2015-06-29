cat("#include <R.h>
     void hello (void) {
     Rprintf(\"hello world!\\n\");
     }",
    file="hello.c")

system("R CMD SHLIB hello.c")

dyn.load(paste0("hello",.Platform$dynlib.ext))

.C("hello",PACKAGE="hello")
