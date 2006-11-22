from distutils.core import setup, Extension 
 
setup(name = "gbfp", 
        version = "0.1", 
        description="gbff parser extension module", 
        author = "Tae-ho Lee", 
        ext_modules=[Extension("gbfp", ["gbfp.python.c", "gbfp.c"])] 
) 
