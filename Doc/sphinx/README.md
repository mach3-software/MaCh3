# Sphinx Documentation

The pyMACh3 module is documented using sphinx. Currently you need to build this documentation yourself if you want it (in future this will be automated). To do this, you will need to install Mach3 with its python extension as described above, then go to the [Doc/sphinx](Doc/sphinx) directory. Then you will need to install sphinx and the necessary extensions which can be done using 

```
pip install -r requirements.txt
```

then you can simply do 

```
make html
```

and the documentation will be built in the build/html directory which you can open with whatever browser you like.
