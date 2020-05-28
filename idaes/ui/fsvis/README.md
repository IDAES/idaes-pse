# IDAES FLOWSHEET VISUALIZER

## For developers
The directory structure follows established conventions for typical medium-sized Flask apps:

````
/fsvis
    /static
        /images
            /icons
    /templates
````

 - The directory `fsvis` contains all python relevant to the `fsvis` module, as well
as source files for running the Flask server used by said module.
 
 - The directory `static` contains all files that are to be used directly by the Flask app
i.e. images, javascript, and any html intended to be served as-is.

 - The directory `templates` holds Flask app templates.
