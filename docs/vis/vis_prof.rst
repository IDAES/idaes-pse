Plotting profile plots from the MEA example
-------------------------------------------

.. warning:: The following has not been tested recently and should be considered a work in progress.

The following examples demonstrate the resize, annotation and saving functionalities.

In the following example, we being by preparing a data frame from our flowsheet variables.

.. code:: python

    # Absorber CO2 Levels
    from pandas import DataFrame
    import os
    tmp = fs.absorb.make_profile(t=0)
    tmp = fs.regen.make_profile(t=0)

    plot_dict = {'z':fs.absorb.profile_1['z'],
                 'y1':fs.absorb.profile_1.y_vap_CO2*101325.0,
                 'y2':fs.absorb.profile_1.P_star_CO2}
    plot_data_frame = DataFrame(data=plot_dict)

We can then plot the data frame we just made, show it, resize it and save it.

.. code:: python

    absorber_co2_plot = Plot.profile(plot_data_frame,
                                     x = 'z',
                                     y = ['y1','y2'],
                                     title = 'Absorber CO2 Levels',
                                     xlab = 'Axial distance from top (m)',
                                     ylab = 'Partial Pressure CO2 (Pa)',
                                     legend = ['Bulk vapor','Equilibrium'])

    absorber_co2_plot.show()
    absorber_co2_plot.save('/home/jovyan/model_contrib/absorber_co2_plot.html')
    assert(os.path.isfile('/home/jovyan/model_contrib/absorber_co2_plot.html'))



.. raw:: html


    <!DOCTYPE html>
    <html lang="en">
        <head>
            <meta charset="utf-8">
            <title>Bokeh Application</title>

    <link rel="stylesheet" href="https://cdn.pydata.org/bokeh/release/bokeh-0.12.9.min.css" type="text/css" />

    <script type="text/javascript" src="https://cdn.pydata.org/bokeh/release/bokeh-0.12.9.min.js"></script>
    <script type="text/javascript">
        Bokeh.set_log_level("info");
    </script>
            <style>
              html {
                width: 100%;
                height: 100%;
              }
              body {
                width: 90%;
                height: 100%;
                margin: auto;
              }
            </style>
        </head>
        <body>

            <div class="bk-root">
                <div class="bk-plotdiv" id="44934bf2-c476-437d-ac6b-508afe9cdbbd"></div>
            </div>

            <script type="text/javascript">
                (function() {
              var fn = function() {
                Bokeh.safely(function() {
                  (function(root) {
                    function embed_document(root) {
                      var docs_json = {"42fca67b-42f1-4540-83b7-8038dd49c5eb":{"roots":{"references":[{"attributes":{"callback":null,"column_names":["x","y"],"data":{"x":{"__ndarray__":"AAAAAAAAAABPOpYZZALKP086lhlkAto/vKswE8uB4z9POpYZZALqP9OuI72RQfA/HnZWQN6B8z9nPYnDKsL2P7EEvEZ3Avo//MvuycNC/T+jyZAmiEEAQEctKmiu4QFA7ZDDqdSBA0CR9Fzr+iEFQGc9icMqwgZADaEiBVFiCECxBLxGdwIKQFZoVYidogtA/MvuycNCDUCgL4gL6uIOQKPJkCaIQRBAdXtdR5sREUBHLSporuERQLNRQFTGsRJAhQMNddmBE0BYtdmV7FEUQCpnprb/IRVA/Bhz1xLyFUDPyj/4JcIWQKJ8DBk5khdAdC7ZOUxiGEBG4KVaXzIZQLEEvEZ3AhpAhLaIZ4rSGkBWaFWInaIbQCgaIqmwchxA/MvuycNCHUDOfbvq1hIeQKAviAvq4h5Ac+FULP2yH0CjyZAmiEEgQNjbmxyUqSBAwTQCrZ0RIUCqjWg9p3khQJTmzs2w4SFAfT81XrpJIkBmmJvuw7EiQE/xAX/NGSNAOEpoD9eBI0Aio86f4OkjQFi12ZXsUSRAQQ5AJva5JEAqZ6a2/yElQBPADEcJiiVA/Bhz1xLyJUDlcdlnHFomQM/KP/glwiZAuSOmiC8qJ0CifAwZOZInQNeOFw9F+idAwOd9n05iKECpQOQvWMooQJKZSsBhMilAe/KwUGuaKUBmSxfhdAIqQE+kfXF+aipAOP3jAYjSKkAhVkqSkTorQFZoVYidoitAP8G7GKcKLEAoGiKpsHIsQBFziDm62ixA/MvuycNCLUDlJFVazaotQM59u+rWEi5At9Yhe+B6LkCgL4gL6uIuQNVBkwH2Si9Avpr5kf+yL0DU+S+RhA0wQEkmY1mJQTBAvVKWIY51MEAyf8npkqkwQKar/LGX3TBAG9gvepwRMUCQBGNCoUUxQKqNaD2neTFAH7qbBaytMUCU5s7NsOExQAgTApa1FTJAfT81XrpJMkA=","dtype":"float64","shape":[91]},"y":{"__ndarray__":"cQovXEbwfkB2hCxogWx/QBq3k0NtsX9AWVnz0UjWf0B9H3jiYul/QI1PM5Mi839AmjXLSR/4f0DERrgmv/p/QORoxSU5/H9AzcvvHSn9f0A3pNwj3f1/QJLUNtp8/n9Angd+Fh3/f0DXDyCgyf9/QMzxDQVFAIBAK9kWBbIAgEAojbZlLgGAQJlkTsG8AYBAw38M6F8CgEAWCGP+GgOAQOC7f5TxA4BA87eUu+cEgEDY134bAgaAQB2xfxFGB4BApRySrrkIgEBVQa79YwqAQG50tAtNDIBAB7CyFH4OgED/BoawARGAQNBhBgbkE4BAxpS0BTMXgEBF1AWt/hqAQDSMgWxZH4BAXbYkIFgkgEBIsxkAEyqAQHkcRrulMIBA8AB9EjA4gEDLP3hx1kCAQPGgNp7CSoBADWwAgiRWgEA11ssPM2OAQKN/C6ItcoBAJN6h3luDgEAWpkzwEJeAQHsvqt6rrYBAVY2an5nHgEA1Ho8fV+WAQJlf/ZRzB4FALBOQKZMugUBvjR0FcluBQDaRke7ojoFAcLXj1ezJgUA9yrfJmQ2CQHGMX/8zW4JAwzsaui60gkAULAEDMxqDQEsEw00nj4NAPUwWKDgVhECSW8AU4q6EQMaNFaQAX4VA0dQ978sohkCgN9TJ/A+HQA0IvabOGIhATETMCxVIiUB42x3AT6OKQOLgSFLBMIxAKoDOJIj3jUBuySEtuv+PQNVgkl5IKZFANM3fpit9kkAWF+rybAGUQIFD79euvJVACh4JCGy2l0B+4V8KD/eZQLZS4C4MiJxA87/rT/5zn0BrW/W3Y2OhQINSfwHmRqNAx6Vk4OhrpUC+b3SE2dqnQOXtnZI1napAbYDZWNW9rUAlNWi1qKSwQAIcZuJTp7JAySVWACnwtEBFXSMQFou3QDc5PdCHiLpA0l2OMpwAvkCW+65Nbg3BQFR/GF5OkMNAzczMzITbxkA=","dtype":"float64","shape":[91]}}},"id":"98c4c9c7-c6c2-401d-b0b1-e08233f15309","type":"ColumnDataSource"},{"attributes":{"dimension":1,"plot":{"id":"64ee0cc1-84cb-4d82-b933-e22da8b399fa","subtype":"Figure","type":"Plot"},"ticker":{"id":"703cdc8c-dbec-4330-9653-d53b4ff2f497","type":"BasicTicker"}},"id":"2e64c291-f3f0-4c18-af6c-435e415d5c13","type":"Grid"},{"attributes":{"line_alpha":{"value":0.1},"line_color":{"value":"#1f77b4"},"x":{"field":"x"},"y":{"field":"y"}},"id":"24f81ed0-7917-4ade-b7f3-841157330ec5","type":"Line"},{"attributes":{},"id":"23a7f411-d56d-43e8-870e-48508bf4b7ce","type":"BasicTickFormatter"},{"attributes":{"line_alpha":{"value":0.1},"line_color":{"value":"#1f77b4"},"x":{"field":"x"},"y":{"field":"y"}},"id":"8ee00533-9f20-402d-a2b1-0b348d98813a","type":"Line"},{"attributes":{"axis_label":"Axial distance from top (m)","formatter":{"id":"54cf35d1-fcce-4715-913f-8e17f7cf7ef4","type":"BasicTickFormatter"},"plot":{"id":"64ee0cc1-84cb-4d82-b933-e22da8b399fa","subtype":"Figure","type":"Plot"},"ticker":{"id":"2dc292fb-ed45-4464-b5a3-c62577f6cf78","type":"BasicTicker"}},"id":"2fff128b-87ec-44d6-b5d9-0354ccb82d4e","type":"LinearAxis"},{"attributes":{"plot":null,"text":"Absorber CO2 Levels"},"id":"e954bc3e-7734-4e68-a76d-5666a4b01b6c","type":"Title"},{"attributes":{},"id":"2dc292fb-ed45-4464-b5a3-c62577f6cf78","type":"BasicTicker"},{"attributes":{"active_drag":"auto","active_inspect":"auto","active_scroll":"auto","active_tap":"auto","tools":[{"id":"d0a8210d-e747-4628-b7bc-461ea42209b4","type":"PanTool"},{"id":"18573688-269b-418a-9ba4-fc63a4bc7c89","type":"BoxZoomTool"},{"id":"8c4571b4-f802-438c-adda-8562942c8e1c","type":"ResetTool"},{"id":"65848d0c-8fcf-4162-91d8-7f320515edbe","type":"SaveTool"}]},"id":"3f44b708-a1e2-4481-b10d-b6996eb88af7","type":"Toolbar"},{"attributes":{},"id":"d0a8210d-e747-4628-b7bc-461ea42209b4","type":"PanTool"},{"attributes":{},"id":"65848d0c-8fcf-4162-91d8-7f320515edbe","type":"SaveTool"},{"attributes":{"line_color":{"value":"green"},"x":{"field":"x"},"y":{"field":"y"}},"id":"4ccea230-59d1-436b-9ba6-e90109e42d9d","type":"Line"},{"attributes":{"source":{"id":"98c4c9c7-c6c2-401d-b0b1-e08233f15309","type":"ColumnDataSource"}},"id":"f1791e55-92cb-49c3-99b5-19cef72585b8","type":"CDSView"},{"attributes":{},"id":"54cf35d1-fcce-4715-913f-8e17f7cf7ef4","type":"BasicTickFormatter"},{"attributes":{"source":{"id":"ca44fe42-60cc-49c3-b388-1db9565e385a","type":"ColumnDataSource"}},"id":"acf0988d-f0d3-49d3-88ca-482b90df7871","type":"CDSView"},{"attributes":{},"id":"2f65a324-84b2-4bf6-949a-34cd82123300","type":"LinearScale"},{"attributes":{},"id":"703cdc8c-dbec-4330-9653-d53b4ff2f497","type":"BasicTicker"},{"attributes":{"label":{"value":"Equilibrium"},"renderers":[{"id":"0ad14117-0a4d-4f26-951d-ef8023bf8e95","type":"GlyphRenderer"}]},"id":"68e529a8-be1d-4d3e-8cd9-f21e7b7a1fcb","type":"LegendItem"},{"attributes":{"line_color":{"value":"blue"},"x":{"field":"x"},"y":{"field":"y"}},"id":"f62183db-489e-4c76-89f5-a1072568fcff","type":"Line"},{"attributes":{"overlay":{"id":"3002211a-5bf2-4297-9225-2ba0482d76e9","type":"BoxAnnotation"}},"id":"18573688-269b-418a-9ba4-fc63a4bc7c89","type":"BoxZoomTool"},{"attributes":{},"id":"8c4571b4-f802-438c-adda-8562942c8e1c","type":"ResetTool"},{"attributes":{"label":{"value":"Bulk vapor"},"renderers":[{"id":"a6245981-f24b-4433-bbe2-6263146e0098","type":"GlyphRenderer"}]},"id":"b523ff34-9b1f-4eb8-a635-ef1923de3b89","type":"LegendItem"},{"attributes":{"plot":{"id":"64ee0cc1-84cb-4d82-b933-e22da8b399fa","subtype":"Figure","type":"Plot"},"ticker":{"id":"2dc292fb-ed45-4464-b5a3-c62577f6cf78","type":"BasicTicker"}},"id":"9fb55311-d186-4c95-88c7-c8e4f5bb26c5","type":"Grid"},{"attributes":{"data_source":{"id":"ca44fe42-60cc-49c3-b388-1db9565e385a","type":"ColumnDataSource"},"glyph":{"id":"f62183db-489e-4c76-89f5-a1072568fcff","type":"Line"},"hover_glyph":null,"muted_glyph":null,"nonselection_glyph":{"id":"8ee00533-9f20-402d-a2b1-0b348d98813a","type":"Line"},"selection_glyph":null,"view":{"id":"acf0988d-f0d3-49d3-88ca-482b90df7871","type":"CDSView"}},"id":"0ad14117-0a4d-4f26-951d-ef8023bf8e95","type":"GlyphRenderer"},{"attributes":{"data_source":{"id":"98c4c9c7-c6c2-401d-b0b1-e08233f15309","type":"ColumnDataSource"},"glyph":{"id":"4ccea230-59d1-436b-9ba6-e90109e42d9d","type":"Line"},"hover_glyph":null,"muted_glyph":null,"nonselection_glyph":{"id":"24f81ed0-7917-4ade-b7f3-841157330ec5","type":"Line"},"selection_glyph":null,"view":{"id":"f1791e55-92cb-49c3-99b5-19cef72585b8","type":"CDSView"}},"id":"a6245981-f24b-4433-bbe2-6263146e0098","type":"GlyphRenderer"},{"attributes":{"below":[{"id":"2fff128b-87ec-44d6-b5d9-0354ccb82d4e","type":"LinearAxis"}],"left":[{"id":"29d33ec9-f0a9-41fb-9377-5779ed729bf8","type":"LinearAxis"}],"renderers":[{"id":"2fff128b-87ec-44d6-b5d9-0354ccb82d4e","type":"LinearAxis"},{"id":"9fb55311-d186-4c95-88c7-c8e4f5bb26c5","type":"Grid"},{"id":"29d33ec9-f0a9-41fb-9377-5779ed729bf8","type":"LinearAxis"},{"id":"2e64c291-f3f0-4c18-af6c-435e415d5c13","type":"Grid"},{"id":"3002211a-5bf2-4297-9225-2ba0482d76e9","type":"BoxAnnotation"},{"id":"ffb375bc-7694-43c3-b98b-c85f8ee4122b","type":"Legend"},{"id":"a6245981-f24b-4433-bbe2-6263146e0098","type":"GlyphRenderer"},{"id":"0ad14117-0a4d-4f26-951d-ef8023bf8e95","type":"GlyphRenderer"}],"title":{"id":"e954bc3e-7734-4e68-a76d-5666a4b01b6c","type":"Title"},"toolbar":{"id":"3f44b708-a1e2-4481-b10d-b6996eb88af7","type":"Toolbar"},"x_range":{"id":"467accd2-8300-44d2-9c43-6aa77ee7b21c","type":"DataRange1d"},"x_scale":{"id":"9237857f-78ed-4ffd-a867-043cbefc51b1","type":"LinearScale"},"y_range":{"id":"3ed75e5b-23b5-4073-b2e3-878d3ba65c49","type":"DataRange1d"},"y_scale":{"id":"2f65a324-84b2-4bf6-949a-34cd82123300","type":"LinearScale"}},"id":"64ee0cc1-84cb-4d82-b933-e22da8b399fa","subtype":"Figure","type":"Plot"},{"attributes":{},"id":"9237857f-78ed-4ffd-a867-043cbefc51b1","type":"LinearScale"},{"attributes":{"items":[{"id":"b523ff34-9b1f-4eb8-a635-ef1923de3b89","type":"LegendItem"},{"id":"68e529a8-be1d-4d3e-8cd9-f21e7b7a1fcb","type":"LegendItem"}],"plot":{"id":"64ee0cc1-84cb-4d82-b933-e22da8b399fa","subtype":"Figure","type":"Plot"}},"id":"ffb375bc-7694-43c3-b98b-c85f8ee4122b","type":"Legend"},{"attributes":{"bottom_units":"screen","fill_alpha":{"value":0.5},"fill_color":{"value":"lightgrey"},"left_units":"screen","level":"overlay","line_alpha":{"value":1.0},"line_color":{"value":"black"},"line_dash":[4,4],"line_width":{"value":2},"plot":null,"render_mode":"css","right_units":"screen","top_units":"screen"},"id":"3002211a-5bf2-4297-9225-2ba0482d76e9","type":"BoxAnnotation"},{"attributes":{"axis_label":"Partial Pressure CO2 (Pa)","formatter":{"id":"23a7f411-d56d-43e8-870e-48508bf4b7ce","type":"BasicTickFormatter"},"plot":{"id":"64ee0cc1-84cb-4d82-b933-e22da8b399fa","subtype":"Figure","type":"Plot"},"ticker":{"id":"703cdc8c-dbec-4330-9653-d53b4ff2f497","type":"BasicTicker"}},"id":"29d33ec9-f0a9-41fb-9377-5779ed729bf8","type":"LinearAxis"},{"attributes":{"callback":null,"column_names":["x","y"],"data":{"x":{"__ndarray__":"AAAAAAAAAABPOpYZZALKP086lhlkAto/vKswE8uB4z9POpYZZALqP9OuI72RQfA/HnZWQN6B8z9nPYnDKsL2P7EEvEZ3Avo//MvuycNC/T+jyZAmiEEAQEctKmiu4QFA7ZDDqdSBA0CR9Fzr+iEFQGc9icMqwgZADaEiBVFiCECxBLxGdwIKQFZoVYidogtA/MvuycNCDUCgL4gL6uIOQKPJkCaIQRBAdXtdR5sREUBHLSporuERQLNRQFTGsRJAhQMNddmBE0BYtdmV7FEUQCpnprb/IRVA/Bhz1xLyFUDPyj/4JcIWQKJ8DBk5khdAdC7ZOUxiGEBG4KVaXzIZQLEEvEZ3AhpAhLaIZ4rSGkBWaFWInaIbQCgaIqmwchxA/MvuycNCHUDOfbvq1hIeQKAviAvq4h5Ac+FULP2yH0CjyZAmiEEgQNjbmxyUqSBAwTQCrZ0RIUCqjWg9p3khQJTmzs2w4SFAfT81XrpJIkBmmJvuw7EiQE/xAX/NGSNAOEpoD9eBI0Aio86f4OkjQFi12ZXsUSRAQQ5AJva5JEAqZ6a2/yElQBPADEcJiiVA/Bhz1xLyJUDlcdlnHFomQM/KP/glwiZAuSOmiC8qJ0CifAwZOZInQNeOFw9F+idAwOd9n05iKECpQOQvWMooQJKZSsBhMilAe/KwUGuaKUBmSxfhdAIqQE+kfXF+aipAOP3jAYjSKkAhVkqSkTorQFZoVYidoitAP8G7GKcKLEAoGiKpsHIsQBFziDm62ixA/MvuycNCLUDlJFVazaotQM59u+rWEi5At9Yhe+B6LkCgL4gL6uIuQNVBkwH2Si9Avpr5kf+yL0DU+S+RhA0wQEkmY1mJQTBAvVKWIY51MEAyf8npkqkwQKar/LGX3TBAG9gvepwRMUCQBGNCoUUxQKqNaD2neTFAH7qbBaytMUCU5s7NsOExQAgTApa1FTJAfT81XrpJMkA=","dtype":"float64","shape":[91]},"y":{"__ndarray__":"wUky3NfPZUC4ztDvKXdxQOYGY3dAxndAWVPI439ifED/APTFNDJ/QF+A217MXoBAVr3X5oTFgEBNQECUZfmAQCmjQP86E4FAXFozxQIggUA+EinAUSaBQG7R6uJvKYFAVYzarPwqgUBdQmTpxCuBQCBjaWgtLIFASa/zzGcsgUAM01abjCyBQAb3xQOoLIFAJ8UnL8AsgUDdEWBX2CyBQLngR03yLIFAZMU9OQ8tgUDpg/b6Ly2BQERsIltVLYFA2ixcIIAtgUDMZKIksS2BQPBKzVnpLYFAaowF0ikugUCGwzLGcy6BQG09gpzILoFASZ947ykvgUCWodmVmS+BQKWdh68ZMIFAGjwEoKwwgUAz0kUyVTGBQBusZZUWMoFAx2+XcPQygUDRxs708jOBQCuMBPEWNYFAA9h76WU2gUBH4HYy5jeBQBOTTxyfOYFAAvtv4Jg7gUBiibcc3T2BQHUbFMh2QIFAFaRXd3JDgUBCrhqa3kaBQEJAssDLSoFAltSG7ExPgUBbCFHsd1SBQB1BYPRlWoFA/TYoXTNhgUBkbElNAWmBQJfNipL1cYFAWF7Djzt8gUAF1WEPBYiBQDMQQjWLlYFAp8kylA+lgUDYozJs3baBQECA+7dLy4FAPYZpRbzigUC+jid1of2BQIqnIad9HIJAEn9OYeY/gkAAeckLh2iCQPt4Dgckl4JA5tnfIJ7MgkCaiLJo9gmDQM97P5VUUINA2uSBhQKhg0BYbMCZfv2DQPesZZt2Z4RAygRlO87ghEDkvZigoWuFQOa+CXBFCoZA9vM9YUK/hkDBgxRVSY2HQOPRsEEjd4hA/i6JbGR/iUBvunUFVKiKQBJvycpL84tACuLLE/ZfjUDZYU1BBeuOQOnO7LMRRpBAg0vf1U8ZkUDsSfJYB+CRQHjvumnUf5JAQS5m+NnLkkDCFhln6nmSQM8YpAi2EpFAUWO8ohPRi0A=","dtype":"float64","shape":[91]}}},"id":"ca44fe42-60cc-49c3-b388-1db9565e385a","type":"ColumnDataSource"},{"attributes":{"callback":null},"id":"467accd2-8300-44d2-9c43-6aa77ee7b21c","type":"DataRange1d"},{"attributes":{"callback":null},"id":"3ed75e5b-23b5-4073-b2e3-878d3ba65c49","type":"DataRange1d"}],"root_ids":["64ee0cc1-84cb-4d82-b933-e22da8b399fa"]},"title":"Bokeh Application","version":"0.12.9"}};
                      var render_items = [{"docid":"42fca67b-42f1-4540-83b7-8038dd49c5eb","elementid":"44934bf2-c476-437d-ac6b-508afe9cdbbd","modelid":"64ee0cc1-84cb-4d82-b933-e22da8b399fa"}];

                      root.Bokeh.embed.embed_items(docs_json, render_items);
                    }

                    if (root.Bokeh !== undefined) {
                      embed_document(root);
                    } else {
                      var attempts = 0;
                      var timer = setInterval(function(root) {
                        if (root.Bokeh !== undefined) {
                          embed_document(root);
                          clearInterval(timer);
                        }
                        attempts++;
                        if (attempts > 100) {
                          console.log("Bokeh: ERROR: Unable to embed document because BokehJS library is missing")
                          clearInterval(timer);
                        }
                      }, 10, root)
                    }
                  })(window);
                });
              };
              if (document.readyState != "loading") fn();
              else document.addEventListener("DOMContentLoaded", fn);
            })();

            </script>
        </body>
    </html>


.. code:: python

    absorber_co2_plot.resize(height=400,width=600)
    absorber_co2_plot.show()
    absorber_co2_plot.save('/home/jovyan/model_contrib/absorber_co2_plot_resized.html')
    assert(os.path.isfile('/home/jovyan/model_contrib/absorber_co2_plot_resized.html'))



.. raw:: html


    <!DOCTYPE html>
    <html lang="en">
        <head>
            <meta charset="utf-8">
            <title>Bokeh Application</title>

    <link rel="stylesheet" href="https://cdn.pydata.org/bokeh/release/bokeh-0.12.9.min.css" type="text/css" />

    <script type="text/javascript" src="https://cdn.pydata.org/bokeh/release/bokeh-0.12.9.min.js"></script>
    <script type="text/javascript">
        Bokeh.set_log_level("info");
    </script>
            <style>
              html {
                width: 100%;
                height: 100%;
              }
              body {
                width: 90%;
                height: 100%;
                margin: auto;
              }
            </style>
        </head>
        <body>

            <div class="bk-root">
                <div class="bk-plotdiv" id="b5654286-c1c3-44a1-85fa-1ea976515cd0"></div>
            </div>

            <script type="text/javascript">
                (function() {
              var fn = function() {
                Bokeh.safely(function() {
                  (function(root) {
                    function embed_document(root) {
                      var docs_json = {"71d47451-1400-4bc3-9da6-5267fc25dbe0":{"roots":{"references":[{"attributes":{"callback":null,"column_names":["x","y"],"data":{"x":{"__ndarray__":"AAAAAAAAAABPOpYZZALKP086lhlkAto/vKswE8uB4z9POpYZZALqP9OuI72RQfA/HnZWQN6B8z9nPYnDKsL2P7EEvEZ3Avo//MvuycNC/T+jyZAmiEEAQEctKmiu4QFA7ZDDqdSBA0CR9Fzr+iEFQGc9icMqwgZADaEiBVFiCECxBLxGdwIKQFZoVYidogtA/MvuycNCDUCgL4gL6uIOQKPJkCaIQRBAdXtdR5sREUBHLSporuERQLNRQFTGsRJAhQMNddmBE0BYtdmV7FEUQCpnprb/IRVA/Bhz1xLyFUDPyj/4JcIWQKJ8DBk5khdAdC7ZOUxiGEBG4KVaXzIZQLEEvEZ3AhpAhLaIZ4rSGkBWaFWInaIbQCgaIqmwchxA/MvuycNCHUDOfbvq1hIeQKAviAvq4h5Ac+FULP2yH0CjyZAmiEEgQNjbmxyUqSBAwTQCrZ0RIUCqjWg9p3khQJTmzs2w4SFAfT81XrpJIkBmmJvuw7EiQE/xAX/NGSNAOEpoD9eBI0Aio86f4OkjQFi12ZXsUSRAQQ5AJva5JEAqZ6a2/yElQBPADEcJiiVA/Bhz1xLyJUDlcdlnHFomQM/KP/glwiZAuSOmiC8qJ0CifAwZOZInQNeOFw9F+idAwOd9n05iKECpQOQvWMooQJKZSsBhMilAe/KwUGuaKUBmSxfhdAIqQE+kfXF+aipAOP3jAYjSKkAhVkqSkTorQFZoVYidoitAP8G7GKcKLEAoGiKpsHIsQBFziDm62ixA/MvuycNCLUDlJFVazaotQM59u+rWEi5At9Yhe+B6LkCgL4gL6uIuQNVBkwH2Si9Avpr5kf+yL0DU+S+RhA0wQEkmY1mJQTBAvVKWIY51MEAyf8npkqkwQKar/LGX3TBAG9gvepwRMUCQBGNCoUUxQKqNaD2neTFAH7qbBaytMUCU5s7NsOExQAgTApa1FTJAfT81XrpJMkA=","dtype":"float64","shape":[91]},"y":{"__ndarray__":"cQovXEbwfkB2hCxogWx/QBq3k0NtsX9AWVnz0UjWf0B9H3jiYul/QI1PM5Mi839AmjXLSR/4f0DERrgmv/p/QORoxSU5/H9AzcvvHSn9f0A3pNwj3f1/QJLUNtp8/n9Angd+Fh3/f0DXDyCgyf9/QMzxDQVFAIBAK9kWBbIAgEAojbZlLgGAQJlkTsG8AYBAw38M6F8CgEAWCGP+GgOAQOC7f5TxA4BA87eUu+cEgEDY134bAgaAQB2xfxFGB4BApRySrrkIgEBVQa79YwqAQG50tAtNDIBAB7CyFH4OgED/BoawARGAQNBhBgbkE4BAxpS0BTMXgEBF1AWt/hqAQDSMgWxZH4BAXbYkIFgkgEBIsxkAEyqAQHkcRrulMIBA8AB9EjA4gEDLP3hx1kCAQPGgNp7CSoBADWwAgiRWgEA11ssPM2OAQKN/C6ItcoBAJN6h3luDgEAWpkzwEJeAQHsvqt6rrYBAVY2an5nHgEA1Ho8fV+WAQJlf/ZRzB4FALBOQKZMugUBvjR0FcluBQDaRke7ojoFAcLXj1ezJgUA9yrfJmQ2CQHGMX/8zW4JAwzsaui60gkAULAEDMxqDQEsEw00nj4NAPUwWKDgVhECSW8AU4q6EQMaNFaQAX4VA0dQ978sohkCgN9TJ/A+HQA0IvabOGIhATETMCxVIiUB42x3AT6OKQOLgSFLBMIxAKoDOJIj3jUBuySEtuv+PQNVgkl5IKZFANM3fpit9kkAWF+rybAGUQIFD79euvJVACh4JCGy2l0B+4V8KD/eZQLZS4C4MiJxA87/rT/5zn0BrW/W3Y2OhQINSfwHmRqNAx6Vk4OhrpUC+b3SE2dqnQOXtnZI1napAbYDZWNW9rUAlNWi1qKSwQAIcZuJTp7JAySVWACnwtEBFXSMQFou3QDc5PdCHiLpA0l2OMpwAvkCW+65Nbg3BQFR/GF5OkMNAzczMzITbxkA=","dtype":"float64","shape":[91]}}},"id":"98c4c9c7-c6c2-401d-b0b1-e08233f15309","type":"ColumnDataSource"},{"attributes":{"dimension":1,"plot":{"id":"64ee0cc1-84cb-4d82-b933-e22da8b399fa","subtype":"Figure","type":"Plot"},"ticker":{"id":"703cdc8c-dbec-4330-9653-d53b4ff2f497","type":"BasicTicker"}},"id":"2e64c291-f3f0-4c18-af6c-435e415d5c13","type":"Grid"},{"attributes":{"line_alpha":{"value":0.1},"line_color":{"value":"#1f77b4"},"x":{"field":"x"},"y":{"field":"y"}},"id":"24f81ed0-7917-4ade-b7f3-841157330ec5","type":"Line"},{"attributes":{},"id":"23a7f411-d56d-43e8-870e-48508bf4b7ce","type":"BasicTickFormatter"},{"attributes":{"line_alpha":{"value":0.1},"line_color":{"value":"#1f77b4"},"x":{"field":"x"},"y":{"field":"y"}},"id":"8ee00533-9f20-402d-a2b1-0b348d98813a","type":"Line"},{"attributes":{"axis_label":"Axial distance from top (m)","formatter":{"id":"54cf35d1-fcce-4715-913f-8e17f7cf7ef4","type":"BasicTickFormatter"},"plot":{"id":"64ee0cc1-84cb-4d82-b933-e22da8b399fa","subtype":"Figure","type":"Plot"},"ticker":{"id":"2dc292fb-ed45-4464-b5a3-c62577f6cf78","type":"BasicTicker"}},"id":"2fff128b-87ec-44d6-b5d9-0354ccb82d4e","type":"LinearAxis"},{"attributes":{"plot":null,"text":"Absorber CO2 Levels"},"id":"e954bc3e-7734-4e68-a76d-5666a4b01b6c","type":"Title"},{"attributes":{},"id":"2dc292fb-ed45-4464-b5a3-c62577f6cf78","type":"BasicTicker"},{"attributes":{"active_drag":"auto","active_inspect":"auto","active_scroll":"auto","active_tap":"auto","tools":[{"id":"d0a8210d-e747-4628-b7bc-461ea42209b4","type":"PanTool"},{"id":"18573688-269b-418a-9ba4-fc63a4bc7c89","type":"BoxZoomTool"},{"id":"8c4571b4-f802-438c-adda-8562942c8e1c","type":"ResetTool"},{"id":"65848d0c-8fcf-4162-91d8-7f320515edbe","type":"SaveTool"}]},"id":"3f44b708-a1e2-4481-b10d-b6996eb88af7","type":"Toolbar"},{"attributes":{},"id":"d0a8210d-e747-4628-b7bc-461ea42209b4","type":"PanTool"},{"attributes":{},"id":"65848d0c-8fcf-4162-91d8-7f320515edbe","type":"SaveTool"},{"attributes":{"line_color":{"value":"green"},"x":{"field":"x"},"y":{"field":"y"}},"id":"4ccea230-59d1-436b-9ba6-e90109e42d9d","type":"Line"},{"attributes":{"source":{"id":"98c4c9c7-c6c2-401d-b0b1-e08233f15309","type":"ColumnDataSource"}},"id":"f1791e55-92cb-49c3-99b5-19cef72585b8","type":"CDSView"},{"attributes":{},"id":"54cf35d1-fcce-4715-913f-8e17f7cf7ef4","type":"BasicTickFormatter"},{"attributes":{"source":{"id":"ca44fe42-60cc-49c3-b388-1db9565e385a","type":"ColumnDataSource"}},"id":"acf0988d-f0d3-49d3-88ca-482b90df7871","type":"CDSView"},{"attributes":{},"id":"2f65a324-84b2-4bf6-949a-34cd82123300","type":"LinearScale"},{"attributes":{},"id":"703cdc8c-dbec-4330-9653-d53b4ff2f497","type":"BasicTicker"},{"attributes":{"label":{"value":"Equilibrium"},"renderers":[{"id":"0ad14117-0a4d-4f26-951d-ef8023bf8e95","type":"GlyphRenderer"}]},"id":"68e529a8-be1d-4d3e-8cd9-f21e7b7a1fcb","type":"LegendItem"},{"attributes":{"line_color":{"value":"blue"},"x":{"field":"x"},"y":{"field":"y"}},"id":"f62183db-489e-4c76-89f5-a1072568fcff","type":"Line"},{"attributes":{"overlay":{"id":"3002211a-5bf2-4297-9225-2ba0482d76e9","type":"BoxAnnotation"}},"id":"18573688-269b-418a-9ba4-fc63a4bc7c89","type":"BoxZoomTool"},{"attributes":{},"id":"8c4571b4-f802-438c-adda-8562942c8e1c","type":"ResetTool"},{"attributes":{"label":{"value":"Bulk vapor"},"renderers":[{"id":"a6245981-f24b-4433-bbe2-6263146e0098","type":"GlyphRenderer"}]},"id":"b523ff34-9b1f-4eb8-a635-ef1923de3b89","type":"LegendItem"},{"attributes":{"plot":{"id":"64ee0cc1-84cb-4d82-b933-e22da8b399fa","subtype":"Figure","type":"Plot"},"ticker":{"id":"2dc292fb-ed45-4464-b5a3-c62577f6cf78","type":"BasicTicker"}},"id":"9fb55311-d186-4c95-88c7-c8e4f5bb26c5","type":"Grid"},{"attributes":{"data_source":{"id":"ca44fe42-60cc-49c3-b388-1db9565e385a","type":"ColumnDataSource"},"glyph":{"id":"f62183db-489e-4c76-89f5-a1072568fcff","type":"Line"},"hover_glyph":null,"muted_glyph":null,"nonselection_glyph":{"id":"8ee00533-9f20-402d-a2b1-0b348d98813a","type":"Line"},"selection_glyph":null,"view":{"id":"acf0988d-f0d3-49d3-88ca-482b90df7871","type":"CDSView"}},"id":"0ad14117-0a4d-4f26-951d-ef8023bf8e95","type":"GlyphRenderer"},{"attributes":{"data_source":{"id":"98c4c9c7-c6c2-401d-b0b1-e08233f15309","type":"ColumnDataSource"},"glyph":{"id":"4ccea230-59d1-436b-9ba6-e90109e42d9d","type":"Line"},"hover_glyph":null,"muted_glyph":null,"nonselection_glyph":{"id":"24f81ed0-7917-4ade-b7f3-841157330ec5","type":"Line"},"selection_glyph":null,"view":{"id":"f1791e55-92cb-49c3-99b5-19cef72585b8","type":"CDSView"}},"id":"a6245981-f24b-4433-bbe2-6263146e0098","type":"GlyphRenderer"},{"attributes":{"below":[{"id":"2fff128b-87ec-44d6-b5d9-0354ccb82d4e","type":"LinearAxis"}],"left":[{"id":"29d33ec9-f0a9-41fb-9377-5779ed729bf8","type":"LinearAxis"}],"plot_height":400,"renderers":[{"id":"2fff128b-87ec-44d6-b5d9-0354ccb82d4e","type":"LinearAxis"},{"id":"9fb55311-d186-4c95-88c7-c8e4f5bb26c5","type":"Grid"},{"id":"29d33ec9-f0a9-41fb-9377-5779ed729bf8","type":"LinearAxis"},{"id":"2e64c291-f3f0-4c18-af6c-435e415d5c13","type":"Grid"},{"id":"3002211a-5bf2-4297-9225-2ba0482d76e9","type":"BoxAnnotation"},{"id":"ffb375bc-7694-43c3-b98b-c85f8ee4122b","type":"Legend"},{"id":"a6245981-f24b-4433-bbe2-6263146e0098","type":"GlyphRenderer"},{"id":"0ad14117-0a4d-4f26-951d-ef8023bf8e95","type":"GlyphRenderer"}],"title":{"id":"e954bc3e-7734-4e68-a76d-5666a4b01b6c","type":"Title"},"toolbar":{"id":"3f44b708-a1e2-4481-b10d-b6996eb88af7","type":"Toolbar"},"x_range":{"id":"467accd2-8300-44d2-9c43-6aa77ee7b21c","type":"DataRange1d"},"x_scale":{"id":"9237857f-78ed-4ffd-a867-043cbefc51b1","type":"LinearScale"},"y_range":{"id":"3ed75e5b-23b5-4073-b2e3-878d3ba65c49","type":"DataRange1d"},"y_scale":{"id":"2f65a324-84b2-4bf6-949a-34cd82123300","type":"LinearScale"}},"id":"64ee0cc1-84cb-4d82-b933-e22da8b399fa","subtype":"Figure","type":"Plot"},{"attributes":{},"id":"9237857f-78ed-4ffd-a867-043cbefc51b1","type":"LinearScale"},{"attributes":{"items":[{"id":"b523ff34-9b1f-4eb8-a635-ef1923de3b89","type":"LegendItem"},{"id":"68e529a8-be1d-4d3e-8cd9-f21e7b7a1fcb","type":"LegendItem"}],"plot":{"id":"64ee0cc1-84cb-4d82-b933-e22da8b399fa","subtype":"Figure","type":"Plot"}},"id":"ffb375bc-7694-43c3-b98b-c85f8ee4122b","type":"Legend"},{"attributes":{"bottom_units":"screen","fill_alpha":{"value":0.5},"fill_color":{"value":"lightgrey"},"left_units":"screen","level":"overlay","line_alpha":{"value":1.0},"line_color":{"value":"black"},"line_dash":[4,4],"line_width":{"value":2},"plot":null,"render_mode":"css","right_units":"screen","top_units":"screen"},"id":"3002211a-5bf2-4297-9225-2ba0482d76e9","type":"BoxAnnotation"},{"attributes":{"axis_label":"Partial Pressure CO2 (Pa)","formatter":{"id":"23a7f411-d56d-43e8-870e-48508bf4b7ce","type":"BasicTickFormatter"},"plot":{"id":"64ee0cc1-84cb-4d82-b933-e22da8b399fa","subtype":"Figure","type":"Plot"},"ticker":{"id":"703cdc8c-dbec-4330-9653-d53b4ff2f497","type":"BasicTicker"}},"id":"29d33ec9-f0a9-41fb-9377-5779ed729bf8","type":"LinearAxis"},{"attributes":{"callback":null,"column_names":["x","y"],"data":{"x":{"__ndarray__":"AAAAAAAAAABPOpYZZALKP086lhlkAto/vKswE8uB4z9POpYZZALqP9OuI72RQfA/HnZWQN6B8z9nPYnDKsL2P7EEvEZ3Avo//MvuycNC/T+jyZAmiEEAQEctKmiu4QFA7ZDDqdSBA0CR9Fzr+iEFQGc9icMqwgZADaEiBVFiCECxBLxGdwIKQFZoVYidogtA/MvuycNCDUCgL4gL6uIOQKPJkCaIQRBAdXtdR5sREUBHLSporuERQLNRQFTGsRJAhQMNddmBE0BYtdmV7FEUQCpnprb/IRVA/Bhz1xLyFUDPyj/4JcIWQKJ8DBk5khdAdC7ZOUxiGEBG4KVaXzIZQLEEvEZ3AhpAhLaIZ4rSGkBWaFWInaIbQCgaIqmwchxA/MvuycNCHUDOfbvq1hIeQKAviAvq4h5Ac+FULP2yH0CjyZAmiEEgQNjbmxyUqSBAwTQCrZ0RIUCqjWg9p3khQJTmzs2w4SFAfT81XrpJIkBmmJvuw7EiQE/xAX/NGSNAOEpoD9eBI0Aio86f4OkjQFi12ZXsUSRAQQ5AJva5JEAqZ6a2/yElQBPADEcJiiVA/Bhz1xLyJUDlcdlnHFomQM/KP/glwiZAuSOmiC8qJ0CifAwZOZInQNeOFw9F+idAwOd9n05iKECpQOQvWMooQJKZSsBhMilAe/KwUGuaKUBmSxfhdAIqQE+kfXF+aipAOP3jAYjSKkAhVkqSkTorQFZoVYidoitAP8G7GKcKLEAoGiKpsHIsQBFziDm62ixA/MvuycNCLUDlJFVazaotQM59u+rWEi5At9Yhe+B6LkCgL4gL6uIuQNVBkwH2Si9Avpr5kf+yL0DU+S+RhA0wQEkmY1mJQTBAvVKWIY51MEAyf8npkqkwQKar/LGX3TBAG9gvepwRMUCQBGNCoUUxQKqNaD2neTFAH7qbBaytMUCU5s7NsOExQAgTApa1FTJAfT81XrpJMkA=","dtype":"float64","shape":[91]},"y":{"__ndarray__":"wUky3NfPZUC4ztDvKXdxQOYGY3dAxndAWVPI439ifED/APTFNDJ/QF+A217MXoBAVr3X5oTFgEBNQECUZfmAQCmjQP86E4FAXFozxQIggUA+EinAUSaBQG7R6uJvKYFAVYzarPwqgUBdQmTpxCuBQCBjaWgtLIFASa/zzGcsgUAM01abjCyBQAb3xQOoLIFAJ8UnL8AsgUDdEWBX2CyBQLngR03yLIFAZMU9OQ8tgUDpg/b6Ly2BQERsIltVLYFA2ixcIIAtgUDMZKIksS2BQPBKzVnpLYFAaowF0ikugUCGwzLGcy6BQG09gpzILoFASZ947ykvgUCWodmVmS+BQKWdh68ZMIFAGjwEoKwwgUAz0kUyVTGBQBusZZUWMoFAx2+XcPQygUDRxs708jOBQCuMBPEWNYFAA9h76WU2gUBH4HYy5jeBQBOTTxyfOYFAAvtv4Jg7gUBiibcc3T2BQHUbFMh2QIFAFaRXd3JDgUBCrhqa3kaBQEJAssDLSoFAltSG7ExPgUBbCFHsd1SBQB1BYPRlWoFA/TYoXTNhgUBkbElNAWmBQJfNipL1cYFAWF7Djzt8gUAF1WEPBYiBQDMQQjWLlYFAp8kylA+lgUDYozJs3baBQECA+7dLy4FAPYZpRbzigUC+jid1of2BQIqnIad9HIJAEn9OYeY/gkAAeckLh2iCQPt4Dgckl4JA5tnfIJ7MgkCaiLJo9gmDQM97P5VUUINA2uSBhQKhg0BYbMCZfv2DQPesZZt2Z4RAygRlO87ghEDkvZigoWuFQOa+CXBFCoZA9vM9YUK/hkDBgxRVSY2HQOPRsEEjd4hA/i6JbGR/iUBvunUFVKiKQBJvycpL84tACuLLE/ZfjUDZYU1BBeuOQOnO7LMRRpBAg0vf1U8ZkUDsSfJYB+CRQHjvumnUf5JAQS5m+NnLkkDCFhln6nmSQM8YpAi2EpFAUWO8ohPRi0A=","dtype":"float64","shape":[91]}}},"id":"ca44fe42-60cc-49c3-b388-1db9565e385a","type":"ColumnDataSource"},{"attributes":{"callback":null},"id":"467accd2-8300-44d2-9c43-6aa77ee7b21c","type":"DataRange1d"},{"attributes":{"callback":null},"id":"3ed75e5b-23b5-4073-b2e3-878d3ba65c49","type":"DataRange1d"}],"root_ids":["64ee0cc1-84cb-4d82-b933-e22da8b399fa"]},"title":"Bokeh Application","version":"0.12.9"}};
                      var render_items = [{"docid":"71d47451-1400-4bc3-9da6-5267fc25dbe0","elementid":"b5654286-c1c3-44a1-85fa-1ea976515cd0","modelid":"64ee0cc1-84cb-4d82-b933-e22da8b399fa"}];

                      root.Bokeh.embed.embed_items(docs_json, render_items);
                    }

                    if (root.Bokeh !== undefined) {
                      embed_document(root);
                    } else {
                      var attempts = 0;
                      var timer = setInterval(function(root) {
                        if (root.Bokeh !== undefined) {
                          embed_document(root);
                          clearInterval(timer);
                        }
                        attempts++;
                        if (attempts > 100) {
                          console.log("Bokeh: ERROR: Unable to embed document because BokehJS library is missing")
                          clearInterval(timer);
                        }
                      }, 10, root)
                    }
                  })(window);
                });
              };
              if (document.readyState != "loading") fn();
              else document.addEventListener("DOMContentLoaded", fn);
            })();

            </script>
        </body>
    </html>

The following demonstrates the annotate functionality by plotting a second plot from the same flowsheet.

.. code:: python

    from IPython.core.display import display,HTML
    stripper_co2_plot = Plot.profile(plot_data_frame,
                                     x = 'z',
                                     y = ['y1','y2'],
                                     title = 'Stripper CO2 Levels',
                                     xlab = 'Axial distance from top (m)',
                                     ylab = 'Partial Pressure CO2 (Pa)',
                                     legend = ['Bulk vapor','Equilibrium'])
    stripper_co2_plot.show()
    stripper_co2_plot.save('/home/jovyan/model_contrib/stripper_co2_plot.html')
    assert(os.path.isfile('/home/jovyan/model_contrib/stripper_co2_plot.html'))



.. raw:: html


    <!DOCTYPE html>
    <html lang="en">
        <head>
            <meta charset="utf-8">
            <title>Bokeh Application</title>

    <link rel="stylesheet" href="https://cdn.pydata.org/bokeh/release/bokeh-0.12.9.min.css" type="text/css" />

    <script type="text/javascript" src="https://cdn.pydata.org/bokeh/release/bokeh-0.12.9.min.js"></script>
    <script type="text/javascript">
        Bokeh.set_log_level("info");
    </script>
            <style>
              html {
                width: 100%;
                height: 100%;
              }
              body {
                width: 90%;
                height: 100%;
                margin: auto;
              }
            </style>
        </head>
        <body>

            <div class="bk-root">
                <div class="bk-plotdiv" id="d46a3a01-a807-4f8e-be99-da6b6aac166a"></div>
            </div>

            <script type="text/javascript">
                (function() {
              var fn = function() {
                Bokeh.safely(function() {
                  (function(root) {
                    function embed_document(root) {
                      var docs_json = {"4a37738b-180e-4323-a3fd-9148852458f0":{"roots":{"references":[{"attributes":{"label":{"value":"Bulk vapor"},"renderers":[{"id":"9512ea20-4566-45a9-aa77-6b3a170237c9","type":"GlyphRenderer"}]},"id":"27dc17d5-2d1c-4db5-a830-edf5679012e2","type":"LegendItem"},{"attributes":{"active_drag":"auto","active_inspect":"auto","active_scroll":"auto","active_tap":"auto","tools":[{"id":"15b976d5-207f-4b91-9e56-580abcf5b752","type":"PanTool"},{"id":"cb5762ed-e9a4-4774-94d2-d66eca83117c","type":"BoxZoomTool"},{"id":"96fd0775-2d63-4d36-9dcc-f8c3f0ca0ca0","type":"ResetTool"},{"id":"61feb7c1-327a-41a9-9487-3b2488b0fcb8","type":"SaveTool"}]},"id":"6fd258fe-34c2-43e3-9581-aab3946e21d0","type":"Toolbar"},{"attributes":{},"id":"96fd0775-2d63-4d36-9dcc-f8c3f0ca0ca0","type":"ResetTool"},{"attributes":{},"id":"61feb7c1-327a-41a9-9487-3b2488b0fcb8","type":"SaveTool"},{"attributes":{"callback":null},"id":"d3d39b15-11c0-4760-9b8a-5eee27de6a51","type":"DataRange1d"},{"attributes":{"label":{"value":"Equilibrium"},"renderers":[{"id":"0753661f-95aa-4388-b6bc-2de50ae27dd9","type":"GlyphRenderer"}]},"id":"32f7a2eb-a43f-442f-814d-8e4817788914","type":"LegendItem"},{"attributes":{"items":[{"id":"27dc17d5-2d1c-4db5-a830-edf5679012e2","type":"LegendItem"},{"id":"32f7a2eb-a43f-442f-814d-8e4817788914","type":"LegendItem"}],"plot":{"id":"8f41f3fb-4ac9-438c-9fc4-526047dfe786","subtype":"Figure","type":"Plot"}},"id":"3bdbbc32-717b-4a00-bf0d-834f224c5daf","type":"Legend"},{"attributes":{"plot":{"id":"8f41f3fb-4ac9-438c-9fc4-526047dfe786","subtype":"Figure","type":"Plot"},"ticker":{"id":"55962ce4-33ed-421b-85e9-c1b7b72fdd2a","type":"BasicTicker"}},"id":"869b8456-0881-490f-837a-6a4145f4b888","type":"Grid"},{"attributes":{},"id":"ea1bf303-86c0-4fe6-bce1-70f28f526db1","type":"BasicTickFormatter"},{"attributes":{"callback":null},"id":"13a72d3f-fabe-419c-978c-f00870aef8b0","type":"DataRange1d"},{"attributes":{},"id":"15b976d5-207f-4b91-9e56-580abcf5b752","type":"PanTool"},{"attributes":{"callback":null,"column_names":["x","y"],"data":{"x":{"__ndarray__":"AAAAAAAAAAB8/stb1DPBP3z+y1vUM9E/u/2xib7N2T98/stb1DPhPyWZhdDigOU/xZh459fN6T9jmGv+zBruPwFMrwrhM/E/0csollta8z+gS6Ih1oD1P3DLG61Qp/c/QEuVOMvN+T8Oyw7ERfT7P2OYa/7MGv4/GozyxKMgAEABTK8K4TMBQOkLbFAeRwJA0csolltaA0C4i+XbmG0EQKBLoiHWgAVAiAtfZxOUBkBwyxutUKcHQBoySkqUughAAvIGkNHNCUDqscPVDuEKQNJxgBtM9AtAuTE9YYkHDUCi8fmmxhoOQImxtuwDLg9AuLg5maAgEECsGBg8P6oQQAFMrwrhMxFA9auNrX+9EUDpC2xQHkcSQN1rSvO80BJA0csolltaE0DFKwc5+uMTQLiL5duYbRRArevDfjf3FECgS6Ih1oAVQPZ+OfB3ChZA6d4XkxaUFkDdPvY1tR0XQNGe1NhTpxdAxf6ye/IwGEC5XpEekboYQKy+b8EvRBlAoB5OZM7NGUCVfiwHbVcaQOqxw9UO4RpA3hGieK1qG0DScYAbTPQbQMXRXr7qfRxAuTE9YYkHHUCtkRsEKJEdQKLx+abGGh5AlVHYSWWkHkCJsbbsAy4fQN7kTbultx9AaSIWL6IgIEBjUoWAcWUgQF2C9NFAqiBAVrJjIxDvIEBR4tJ03zMhQEsSQsaueCFARUKxF369IUA+ciBpTQIiQOkLbFAeRyJA4zvboe2LIkDda0rzvNAiQNebuUSMFSNA0csolltaI0DL+5fnKp8jQMUrBzn64yNAv1t2iskoJEC4i+XbmG0kQGMlMcNpsiRAXVWgFDn3JEBXhQ9mCDwlQFG1frfXgCVAS+XtCKfFJUBFFV1adgomQD9FzKtFTyZAOXU7/RSUJkAypapO5NgmQN0+9jW1HSdA2G5lh4RiJ0DRntTYU6cnQMvOQyoj7CdAxf6ye/IwKEA=","dtype":"float64","shape":[91]},"y":{"__ndarray__":"cLAjBGuz7kC1NIFlJJDtQEjUB4lS6+xAc//HbUSP7EAO4krLK1zsQKo2Z9/gP+xA0PWmFTUw7ECJwmO+fifsQGDz9pebIuxA7dVvdNAf7EBdhEaXKB7sQKpwg6ocHexACd2ni2Ec7EDzVF8nzRvsQMOou4BHG+xAbZx4dsIa7EBAluwpNRrsQI09532ZGexAr/spq+oY7EC+/9dyJBjsQM+7DqZCF+xAEiVb20AW7EBACkw8GhXsQAsytlnJE+xAPuWOHkgS7ECpi4uTjxDsQDN24M+XDuxAACWn0VcM7ECMVNBWxQnsQB1S0rHUBuxA81B4mXgD7EDrv0/zof/rQHKcjIY/++tA0d13/T3260ATr5Q+h/DrQBGQYlkC6utAxDzpI5Pi60Akq+7fGdrrQM2CINty0OtAvmKPC3bF60Agyomq9rjrQPVuQpfCqutAbZLi0aKa60Cu127yWYjrQPn7VKakc+tA8JiyWTlc60AE1cQvyEHrQLOuxSn7I+tAA08PiXYC60A++VJ62dzqQFhCiWy+supA/UInA7+D6kBjjB1Jck/qQC/GRj1xFepAfgwUA1nV6UAFRzD9zY7pQNlWkjR/QelAKD+L5int6EB161UHnZHoQD+1s+K6LuhAHXxt8YHE50AH2F7cBVPnQLKAQed12uZAx7ZjRhtb5kAkcHZ8V9XlQPrW8IWhSeVAdG8H/oG45EAwe1iCjiLkQPowNxBiiONAysCSGaLq4kCqj1/S6kniQJlJ9mDUpuFA3u64B+0B4UA+4sJNtlvgQG2fK+lFad9AiHudUisa3kARNPyowMrcQD0MBsCBe9tADBB5ltos2kCQO1C7BN/YQO6MjQYfktdAodm35yhG1kCwKWzWAvvUQP5Kv2FusNNAyHfgxw1m0kCMn2gpYxvRQCgHQ+CSn89A5Q0+nBwFzUB8Fv4UpWXKQJ3w0Xtsv8dAHh84wg8RxUA=","dtype":"float64","shape":[91]}}},"id":"69ab0692-9132-4bb6-af79-570a5c3cdaa6","type":"ColumnDataSource"},{"attributes":{"line_color":{"value":"blue"},"x":{"field":"x"},"y":{"field":"y"}},"id":"7bcce2c0-0680-4927-b068-e68d78fd7202","type":"Line"},{"attributes":{},"id":"bcc2034d-131e-4153-b7b9-d98f10a07037","type":"BasicTickFormatter"},{"attributes":{"line_alpha":{"value":0.1},"line_color":{"value":"#1f77b4"},"x":{"field":"x"},"y":{"field":"y"}},"id":"03c8c8d3-2a0b-4ff5-ba11-419aa59d6ed5","type":"Line"},{"attributes":{"line_color":{"value":"green"},"x":{"field":"x"},"y":{"field":"y"}},"id":"b32732a1-dd3a-48f6-9bd6-5c193983cce3","type":"Line"},{"attributes":{"dimension":1,"plot":{"id":"8f41f3fb-4ac9-438c-9fc4-526047dfe786","subtype":"Figure","type":"Plot"},"ticker":{"id":"702fd5c2-b7c4-46f4-ac99-a654b7a61a1e","type":"BasicTicker"}},"id":"97ec3a3c-01c7-4912-8e62-98284d63dae3","type":"Grid"},{"attributes":{"source":{"id":"c4414caf-b280-422b-8483-1807c1d68dfe","type":"ColumnDataSource"}},"id":"795116a0-8abf-46eb-99bf-a882e1d1b744","type":"CDSView"},{"attributes":{},"id":"cd033eeb-4342-4644-8a81-6e4a57836b62","type":"LinearScale"},{"attributes":{"line_alpha":{"value":0.1},"line_color":{"value":"#1f77b4"},"x":{"field":"x"},"y":{"field":"y"}},"id":"f02a2f3f-3c10-4e60-b23e-3cddc2528ddf","type":"Line"},{"attributes":{"bottom_units":"screen","fill_alpha":{"value":0.5},"fill_color":{"value":"lightgrey"},"left_units":"screen","level":"overlay","line_alpha":{"value":1.0},"line_color":{"value":"black"},"line_dash":[4,4],"line_width":{"value":2},"plot":null,"render_mode":"css","right_units":"screen","top_units":"screen"},"id":"9be237e4-696b-402b-a816-cca962343a6f","type":"BoxAnnotation"},{"attributes":{"plot":null,"text":"Stripper CO2 Levels"},"id":"9135a977-2a41-4547-92f8-3a2959da633a","type":"Title"},{"attributes":{"data_source":{"id":"69ab0692-9132-4bb6-af79-570a5c3cdaa6","type":"ColumnDataSource"},"glyph":{"id":"b32732a1-dd3a-48f6-9bd6-5c193983cce3","type":"Line"},"hover_glyph":null,"muted_glyph":null,"nonselection_glyph":{"id":"03c8c8d3-2a0b-4ff5-ba11-419aa59d6ed5","type":"Line"},"selection_glyph":null,"view":{"id":"ef9ffe0b-4431-4440-a631-3d83902fcfec","type":"CDSView"}},"id":"9512ea20-4566-45a9-aa77-6b3a170237c9","type":"GlyphRenderer"},{"attributes":{"source":{"id":"69ab0692-9132-4bb6-af79-570a5c3cdaa6","type":"ColumnDataSource"}},"id":"ef9ffe0b-4431-4440-a631-3d83902fcfec","type":"CDSView"},{"attributes":{"below":[{"id":"18ba8f4c-a5ea-477a-931e-d58496934361","type":"LinearAxis"}],"left":[{"id":"40860ef2-c2cc-40ed-a2c5-15b19d7116d5","type":"LinearAxis"}],"renderers":[{"id":"18ba8f4c-a5ea-477a-931e-d58496934361","type":"LinearAxis"},{"id":"869b8456-0881-490f-837a-6a4145f4b888","type":"Grid"},{"id":"40860ef2-c2cc-40ed-a2c5-15b19d7116d5","type":"LinearAxis"},{"id":"97ec3a3c-01c7-4912-8e62-98284d63dae3","type":"Grid"},{"id":"9be237e4-696b-402b-a816-cca962343a6f","type":"BoxAnnotation"},{"id":"3bdbbc32-717b-4a00-bf0d-834f224c5daf","type":"Legend"},{"id":"9512ea20-4566-45a9-aa77-6b3a170237c9","type":"GlyphRenderer"},{"id":"0753661f-95aa-4388-b6bc-2de50ae27dd9","type":"GlyphRenderer"}],"title":{"id":"9135a977-2a41-4547-92f8-3a2959da633a","type":"Title"},"toolbar":{"id":"6fd258fe-34c2-43e3-9581-aab3946e21d0","type":"Toolbar"},"x_range":{"id":"d3d39b15-11c0-4760-9b8a-5eee27de6a51","type":"DataRange1d"},"x_scale":{"id":"7090e1e7-53a1-4532-b87c-5e8e94868d5c","type":"LinearScale"},"y_range":{"id":"13a72d3f-fabe-419c-978c-f00870aef8b0","type":"DataRange1d"},"y_scale":{"id":"cd033eeb-4342-4644-8a81-6e4a57836b62","type":"LinearScale"}},"id":"8f41f3fb-4ac9-438c-9fc4-526047dfe786","subtype":"Figure","type":"Plot"},{"attributes":{},"id":"7090e1e7-53a1-4532-b87c-5e8e94868d5c","type":"LinearScale"},{"attributes":{},"id":"55962ce4-33ed-421b-85e9-c1b7b72fdd2a","type":"BasicTicker"},{"attributes":{},"id":"702fd5c2-b7c4-46f4-ac99-a654b7a61a1e","type":"BasicTicker"},{"attributes":{"overlay":{"id":"9be237e4-696b-402b-a816-cca962343a6f","type":"BoxAnnotation"}},"id":"cb5762ed-e9a4-4774-94d2-d66eca83117c","type":"BoxZoomTool"},{"attributes":{"axis_label":"Partial Pressure CO2 (Pa)","formatter":{"id":"ea1bf303-86c0-4fe6-bce1-70f28f526db1","type":"BasicTickFormatter"},"plot":{"id":"8f41f3fb-4ac9-438c-9fc4-526047dfe786","subtype":"Figure","type":"Plot"},"ticker":{"id":"702fd5c2-b7c4-46f4-ac99-a654b7a61a1e","type":"BasicTicker"}},"id":"40860ef2-c2cc-40ed-a2c5-15b19d7116d5","type":"LinearAxis"},{"attributes":{"axis_label":"Axial distance from top (m)","formatter":{"id":"bcc2034d-131e-4153-b7b9-d98f10a07037","type":"BasicTickFormatter"},"plot":{"id":"8f41f3fb-4ac9-438c-9fc4-526047dfe786","subtype":"Figure","type":"Plot"},"ticker":{"id":"55962ce4-33ed-421b-85e9-c1b7b72fdd2a","type":"BasicTicker"}},"id":"18ba8f4c-a5ea-477a-931e-d58496934361","type":"LinearAxis"},{"attributes":{"data_source":{"id":"c4414caf-b280-422b-8483-1807c1d68dfe","type":"ColumnDataSource"},"glyph":{"id":"7bcce2c0-0680-4927-b068-e68d78fd7202","type":"Line"},"hover_glyph":null,"muted_glyph":null,"nonselection_glyph":{"id":"f02a2f3f-3c10-4e60-b23e-3cddc2528ddf","type":"Line"},"selection_glyph":null,"view":{"id":"795116a0-8abf-46eb-99bf-a882e1d1b744","type":"CDSView"}},"id":"0753661f-95aa-4388-b6bc-2de50ae27dd9","type":"GlyphRenderer"},{"attributes":{"callback":null,"column_names":["x","y"],"data":{"x":{"__ndarray__":"AAAAAAAAAAB8/stb1DPBP3z+y1vUM9E/u/2xib7N2T98/stb1DPhPyWZhdDigOU/xZh459fN6T9jmGv+zBruPwFMrwrhM/E/0csollta8z+gS6Ih1oD1P3DLG61Qp/c/QEuVOMvN+T8Oyw7ERfT7P2OYa/7MGv4/GozyxKMgAEABTK8K4TMBQOkLbFAeRwJA0csolltaA0C4i+XbmG0EQKBLoiHWgAVAiAtfZxOUBkBwyxutUKcHQBoySkqUughAAvIGkNHNCUDqscPVDuEKQNJxgBtM9AtAuTE9YYkHDUCi8fmmxhoOQImxtuwDLg9AuLg5maAgEECsGBg8P6oQQAFMrwrhMxFA9auNrX+9EUDpC2xQHkcSQN1rSvO80BJA0csolltaE0DFKwc5+uMTQLiL5duYbRRArevDfjf3FECgS6Ih1oAVQPZ+OfB3ChZA6d4XkxaUFkDdPvY1tR0XQNGe1NhTpxdAxf6ye/IwGEC5XpEekboYQKy+b8EvRBlAoB5OZM7NGUCVfiwHbVcaQOqxw9UO4RpA3hGieK1qG0DScYAbTPQbQMXRXr7qfRxAuTE9YYkHHUCtkRsEKJEdQKLx+abGGh5AlVHYSWWkHkCJsbbsAy4fQN7kTbultx9AaSIWL6IgIEBjUoWAcWUgQF2C9NFAqiBAVrJjIxDvIEBR4tJ03zMhQEsSQsaueCFARUKxF369IUA+ciBpTQIiQOkLbFAeRyJA4zvboe2LIkDda0rzvNAiQNebuUSMFSNA0csolltaI0DL+5fnKp8jQMUrBzn64yNAv1t2iskoJEC4i+XbmG0kQGMlMcNpsiRAXVWgFDn3JEBXhQ9mCDwlQFG1frfXgCVAS+XtCKfFJUBFFV1adgomQD9FzKtFTyZAOXU7/RSUJkAypapO5NgmQN0+9jW1HSdA2G5lh4RiJ0DRntTYU6cnQMvOQyoj7CdAxf6ye/IwKEA=","dtype":"float64","shape":[91]},"y":{"__ndarray__":"tPRumMz54UB+g/DBIQnmQD04nW2/oehAfJJ3Qd8v6kAsqlQslxXrQBJmdERdl+tA8SDmX8zf60AS7nfW9AfsQPMoAbogHuxAi8BboVQq7EA3H70LBjHsQL5c6KOsNOxAPqdR8aQ27ECvdCDhrjfsQNJoFvIzOOxAMnzX92047ED6o6CAfDjsQNjuppdwOOxAJJpzQFM47ED8z9oFKTjsQH64Wu3zN+xAf+qmh7Q37EBrjQyEajfsQD8WO/oUN+xAq6RopLI27ED+cULTQTbsQCB7e4nANexAIpBKcyw17EBRu33fgjTsQHK8/bTAM+xAkd9hZeIy7ECMj8fc4zHsQD+FM2DAMOxABS7btXIv7EB+thWn9C3sQD1/oi0/LOxAWXPHPEoq7EBh0UqbDCjsQF9+0Ld7JexAWR7ldosi7EAAiOr5LR/sQOLdeS1TG+xA2IAfTekW7ECzhbZC2xHsQCABIT0RDOxAM+8I/28F7EC7MEBp2P3rQBIyt/gm9etA0My+ODPr60As3jsrz9/rQDcBWBDG0utAY5s9Tt7D60BtbZS+1LLrQHfGfyRfn+tAwDz13iqJ60A730Wr3G/rQGFgp6sQU+tAzpvUyloy60BlIECXRw3rQMntsAdc4+pATQBvgR606kCZ8FlyDn/qQO98N+2vQ+pAh2w9Vo0B6kB/8qjfO7jpQDRdETFgZ+lAOogy6bIO6UDVOZ2fBK7oQHAabrc9RehA677FTG7U50DrtEyxuFvnQGcbzLxg2+ZAmzc1/MRT5kCrGyJgW8XlQIMeBvWsMOVA6W3P/FCW5EDhEJnC5/bjQMUNyBsRU+NAOjBqt3mr4kBjGd8dugDiQNvbfJRnU+FA/F1VEguk4ECfgoiePebfQC1AtmkYgt5AHCdd31Mc3UDY2LMOc7XbQH1cxKK9TdpA67MWi1Dl2ECfyd3zuXvXQEkx3Ee0D9ZAaeGG3hae1EA=","dtype":"float64","shape":[91]}}},"id":"c4414caf-b280-422b-8483-1807c1d68dfe","type":"ColumnDataSource"}],"root_ids":["8f41f3fb-4ac9-438c-9fc4-526047dfe786"]},"title":"Bokeh Application","version":"0.12.9"}};
                      var render_items = [{"docid":"4a37738b-180e-4323-a3fd-9148852458f0","elementid":"d46a3a01-a807-4f8e-be99-da6b6aac166a","modelid":"8f41f3fb-4ac9-438c-9fc4-526047dfe786"}];

                      root.Bokeh.embed.embed_items(docs_json, render_items);
                    }

                    if (root.Bokeh !== undefined) {
                      embed_document(root);
                    } else {
                      var attempts = 0;
                      var timer = setInterval(function(root) {
                        if (root.Bokeh !== undefined) {
                          embed_document(root);
                          clearInterval(timer);
                        }
                        attempts++;
                        if (attempts > 100) {
                          console.log("Bokeh: ERROR: Unable to embed document because BokehJS library is missing")
                          clearInterval(timer);
                        }
                      }, 10, root)
                    }
                  })(window);
                });
              };
              if (document.readyState != "loading") fn();
              else document.addEventListener("DOMContentLoaded", fn);
            })();

            </script>
        </body>
    </html>


We can then annotate the "Reboiler vapor" point as shown below:

.. code:: python

    stripper_co2_plot.annotate(rloc,rco2p,'Reboiler vapor')
    stripper_co2_plot.show()
    stripper_co2_plot.save('/home/jovyan/model_contrib/stripper_co2_plot_annotated.html')



.. raw:: html


    <!DOCTYPE html>
    <html lang="en">
        <head>
            <meta charset="utf-8">
            <title>Bokeh Application</title>

    <link rel="stylesheet" href="https://cdn.pydata.org/bokeh/release/bokeh-0.12.9.min.css" type="text/css" />

    <script type="text/javascript" src="https://cdn.pydata.org/bokeh/release/bokeh-0.12.9.min.js"></script>
    <script type="text/javascript">
        Bokeh.set_log_level("info");
    </script>
            <style>
              html {
                width: 100%;
                height: 100%;
              }
              body {
                width: 90%;
                height: 100%;
                margin: auto;
              }
            </style>
        </head>
        <body>

            <div class="bk-root">
                <div class="bk-plotdiv" id="e3811af1-a8f9-4042-b961-7400b192ed89"></div>
            </div>

            <script type="text/javascript">
                (function() {
              var fn = function() {
                Bokeh.safely(function() {
                  (function(root) {
                    function embed_document(root) {
                      var docs_json = {"80139986-1908-4271-886f-f4e38b644bfb":{"roots":{"references":[{"attributes":{"callback":null,"column_names":["x","y"],"data":{"x":{"__ndarray__":"AAAAAAAAAAB8/stb1DPBP3z+y1vUM9E/u/2xib7N2T98/stb1DPhPyWZhdDigOU/xZh459fN6T9jmGv+zBruPwFMrwrhM/E/0csollta8z+gS6Ih1oD1P3DLG61Qp/c/QEuVOMvN+T8Oyw7ERfT7P2OYa/7MGv4/GozyxKMgAEABTK8K4TMBQOkLbFAeRwJA0csolltaA0C4i+XbmG0EQKBLoiHWgAVAiAtfZxOUBkBwyxutUKcHQBoySkqUughAAvIGkNHNCUDqscPVDuEKQNJxgBtM9AtAuTE9YYkHDUCi8fmmxhoOQImxtuwDLg9AuLg5maAgEECsGBg8P6oQQAFMrwrhMxFA9auNrX+9EUDpC2xQHkcSQN1rSvO80BJA0csolltaE0DFKwc5+uMTQLiL5duYbRRArevDfjf3FECgS6Ih1oAVQPZ+OfB3ChZA6d4XkxaUFkDdPvY1tR0XQNGe1NhTpxdAxf6ye/IwGEC5XpEekboYQKy+b8EvRBlAoB5OZM7NGUCVfiwHbVcaQOqxw9UO4RpA3hGieK1qG0DScYAbTPQbQMXRXr7qfRxAuTE9YYkHHUCtkRsEKJEdQKLx+abGGh5AlVHYSWWkHkCJsbbsAy4fQN7kTbultx9AaSIWL6IgIEBjUoWAcWUgQF2C9NFAqiBAVrJjIxDvIEBR4tJ03zMhQEsSQsaueCFARUKxF369IUA+ciBpTQIiQOkLbFAeRyJA4zvboe2LIkDda0rzvNAiQNebuUSMFSNA0csolltaI0DL+5fnKp8jQMUrBzn64yNAv1t2iskoJEC4i+XbmG0kQGMlMcNpsiRAXVWgFDn3JEBXhQ9mCDwlQFG1frfXgCVAS+XtCKfFJUBFFV1adgomQD9FzKtFTyZAOXU7/RSUJkAypapO5NgmQN0+9jW1HSdA2G5lh4RiJ0DRntTYU6cnQMvOQyoj7CdAxf6ye/IwKEA=","dtype":"float64","shape":[91]},"y":{"__ndarray__":"cLAjBGuz7kC1NIFlJJDtQEjUB4lS6+xAc//HbUSP7EAO4krLK1zsQKo2Z9/gP+xA0PWmFTUw7ECJwmO+fifsQGDz9pebIuxA7dVvdNAf7EBdhEaXKB7sQKpwg6ocHexACd2ni2Ec7EDzVF8nzRvsQMOou4BHG+xAbZx4dsIa7EBAluwpNRrsQI09532ZGexAr/spq+oY7EC+/9dyJBjsQM+7DqZCF+xAEiVb20AW7EBACkw8GhXsQAsytlnJE+xAPuWOHkgS7ECpi4uTjxDsQDN24M+XDuxAACWn0VcM7ECMVNBWxQnsQB1S0rHUBuxA81B4mXgD7EDrv0/zof/rQHKcjIY/++tA0d13/T3260ATr5Q+h/DrQBGQYlkC6utAxDzpI5Pi60Akq+7fGdrrQM2CINty0OtAvmKPC3bF60Agyomq9rjrQPVuQpfCqutAbZLi0aKa60Cu127yWYjrQPn7VKakc+tA8JiyWTlc60AE1cQvyEHrQLOuxSn7I+tAA08PiXYC60A++VJ62dzqQFhCiWy+supA/UInA7+D6kBjjB1Jck/qQC/GRj1xFepAfgwUA1nV6UAFRzD9zY7pQNlWkjR/QelAKD+L5int6EB161UHnZHoQD+1s+K6LuhAHXxt8YHE50AH2F7cBVPnQLKAQed12uZAx7ZjRhtb5kAkcHZ8V9XlQPrW8IWhSeVAdG8H/oG45EAwe1iCjiLkQPowNxBiiONAysCSGaLq4kCqj1/S6kniQJlJ9mDUpuFA3u64B+0B4UA+4sJNtlvgQG2fK+lFad9AiHudUisa3kARNPyowMrcQD0MBsCBe9tADBB5ltos2kCQO1C7BN/YQO6MjQYfktdAodm35yhG1kCwKWzWAvvUQP5Kv2FusNNAyHfgxw1m0kCMn2gpYxvRQCgHQ+CSn89A5Q0+nBwFzUB8Fv4UpWXKQJ3w0Xtsv8dAHh84wg8RxUA=","dtype":"float64","shape":[91]}}},"id":"69ab0692-9132-4bb6-af79-570a5c3cdaa6","type":"ColumnDataSource"},{"attributes":{"active_drag":"auto","active_inspect":"auto","active_scroll":"auto","active_tap":"auto","tools":[{"id":"15b976d5-207f-4b91-9e56-580abcf5b752","type":"PanTool"},{"id":"cb5762ed-e9a4-4774-94d2-d66eca83117c","type":"BoxZoomTool"},{"id":"96fd0775-2d63-4d36-9dcc-f8c3f0ca0ca0","type":"ResetTool"},{"id":"61feb7c1-327a-41a9-9487-3b2488b0fcb8","type":"SaveTool"}]},"id":"6fd258fe-34c2-43e3-9581-aab3946e21d0","type":"Toolbar"},{"attributes":{"callback":null,"data":{}},"id":"14a6fc3b-7f92-488c-81d4-fe802b1b1b15","type":"ColumnDataSource"},{"attributes":{},"id":"61feb7c1-327a-41a9-9487-3b2488b0fcb8","type":"SaveTool"},{"attributes":{"callback":null},"id":"d3d39b15-11c0-4760-9b8a-5eee27de6a51","type":"DataRange1d"},{"attributes":{"data_source":{"id":"14a6fc3b-7f92-488c-81d4-fe802b1b1b15","type":"ColumnDataSource"},"glyph":{"id":"944d5d75-3a77-4f93-9f55-51c2201afcee","type":"Circle"},"hover_glyph":null,"muted_glyph":null,"nonselection_glyph":{"id":"de2ec76b-81fd-4428-a7ba-cc1c348a021b","type":"Circle"},"selection_glyph":null,"view":{"id":"5cb07060-2988-4b92-b7a5-a3595b59e9e7","type":"CDSView"}},"id":"6d8ec9f4-f4f2-4590-94e2-4170ef2cd166","type":"GlyphRenderer"},{"attributes":{"label":{"value":"Equilibrium"},"renderers":[{"id":"0753661f-95aa-4388-b6bc-2de50ae27dd9","type":"GlyphRenderer"}]},"id":"32f7a2eb-a43f-442f-814d-8e4817788914","type":"LegendItem"},{"attributes":{"bottom_units":"screen","fill_alpha":{"value":0.5},"fill_color":{"value":"lightgrey"},"left_units":"screen","level":"overlay","line_alpha":{"value":1.0},"line_color":{"value":"black"},"line_dash":[4,4],"line_width":{"value":2},"plot":null,"render_mode":"css","right_units":"screen","top_units":"screen"},"id":"9be237e4-696b-402b-a816-cca962343a6f","type":"BoxAnnotation"},{"attributes":{"plot":{"id":"8f41f3fb-4ac9-438c-9fc4-526047dfe786","subtype":"Figure","type":"Plot"},"ticker":{"id":"55962ce4-33ed-421b-85e9-c1b7b72fdd2a","type":"BasicTicker"}},"id":"869b8456-0881-490f-837a-6a4145f4b888","type":"Grid"},{"attributes":{},"id":"ea1bf303-86c0-4fe6-bce1-70f28f526db1","type":"BasicTickFormatter"},{"attributes":{"callback":null},"id":"13a72d3f-fabe-419c-978c-f00870aef8b0","type":"DataRange1d"},{"attributes":{},"id":"15b976d5-207f-4b91-9e56-580abcf5b752","type":"PanTool"},{"attributes":{"line_color":{"value":"blue"},"x":{"field":"x"},"y":{"field":"y"}},"id":"7bcce2c0-0680-4927-b068-e68d78fd7202","type":"Line"},{"attributes":{},"id":"bcc2034d-131e-4153-b7b9-d98f10a07037","type":"BasicTickFormatter"},{"attributes":{"dimension":1,"plot":{"id":"8f41f3fb-4ac9-438c-9fc4-526047dfe786","subtype":"Figure","type":"Plot"},"ticker":{"id":"702fd5c2-b7c4-46f4-ac99-a654b7a61a1e","type":"BasicTicker"}},"id":"97ec3a3c-01c7-4912-8e62-98284d63dae3","type":"Grid"},{"attributes":{"fill_alpha":{"value":0.1},"fill_color":{"value":"#1f77b4"},"line_alpha":{"value":0.1},"line_color":{"value":"#1f77b4"},"x":{"value":12.0956},"y":{"value":10827.3815090278}},"id":"de2ec76b-81fd-4428-a7ba-cc1c348a021b","type":"Circle"},{"attributes":{"line_alpha":{"value":0.1},"line_color":{"value":"#1f77b4"},"x":{"field":"x"},"y":{"field":"y"}},"id":"03c8c8d3-2a0b-4ff5-ba11-419aa59d6ed5","type":"Line"},{"attributes":{"fill_color":{"value":"red"},"line_color":{"value":"red"},"x":{"value":12.0956},"y":{"value":10827.3815090278}},"id":"944d5d75-3a77-4f93-9f55-51c2201afcee","type":"Circle"},{"attributes":{"source":{"id":"c4414caf-b280-422b-8483-1807c1d68dfe","type":"ColumnDataSource"}},"id":"795116a0-8abf-46eb-99bf-a882e1d1b744","type":"CDSView"},{"attributes":{"level":"glyph","plot":{"id":"8f41f3fb-4ac9-438c-9fc4-526047dfe786","subtype":"Figure","type":"Plot"},"source":{"id":"8346089e-d296-4764-91ec-ff4ef05a917b","type":"ColumnDataSource"},"text":{"field":"labels"},"text_font_size":{"value":"9pt"},"x":{"field":"x"},"y":{"field":"y"}},"id":"4088aadd-2cbd-466f-a523-a9c87a1566b0","type":"LabelSet"},{"attributes":{"line_color":{"value":"green"},"x":{"field":"x"},"y":{"field":"y"}},"id":"b32732a1-dd3a-48f6-9bd6-5c193983cce3","type":"Line"},{"attributes":{},"id":"cd033eeb-4342-4644-8a81-6e4a57836b62","type":"LinearScale"},{"attributes":{"line_alpha":{"value":0.1},"line_color":{"value":"#1f77b4"},"x":{"field":"x"},"y":{"field":"y"}},"id":"f02a2f3f-3c10-4e60-b23e-3cddc2528ddf","type":"Line"},{"attributes":{"items":[{"id":"27dc17d5-2d1c-4db5-a830-edf5679012e2","type":"LegendItem"},{"id":"32f7a2eb-a43f-442f-814d-8e4817788914","type":"LegendItem"}],"plot":{"id":"8f41f3fb-4ac9-438c-9fc4-526047dfe786","subtype":"Figure","type":"Plot"}},"id":"3bdbbc32-717b-4a00-bf0d-834f224c5daf","type":"Legend"},{"attributes":{"plot":null,"text":"Stripper CO2 Levels"},"id":"9135a977-2a41-4547-92f8-3a2959da633a","type":"Title"},{"attributes":{"data_source":{"id":"69ab0692-9132-4bb6-af79-570a5c3cdaa6","type":"ColumnDataSource"},"glyph":{"id":"b32732a1-dd3a-48f6-9bd6-5c193983cce3","type":"Line"},"hover_glyph":null,"muted_glyph":null,"nonselection_glyph":{"id":"03c8c8d3-2a0b-4ff5-ba11-419aa59d6ed5","type":"Line"},"selection_glyph":null,"view":{"id":"ef9ffe0b-4431-4440-a631-3d83902fcfec","type":"CDSView"}},"id":"9512ea20-4566-45a9-aa77-6b3a170237c9","type":"GlyphRenderer"},{"attributes":{"source":{"id":"69ab0692-9132-4bb6-af79-570a5c3cdaa6","type":"ColumnDataSource"}},"id":"ef9ffe0b-4431-4440-a631-3d83902fcfec","type":"CDSView"},{"attributes":{"below":[{"id":"18ba8f4c-a5ea-477a-931e-d58496934361","type":"LinearAxis"}],"left":[{"id":"40860ef2-c2cc-40ed-a2c5-15b19d7116d5","type":"LinearAxis"}],"renderers":[{"id":"18ba8f4c-a5ea-477a-931e-d58496934361","type":"LinearAxis"},{"id":"869b8456-0881-490f-837a-6a4145f4b888","type":"Grid"},{"id":"40860ef2-c2cc-40ed-a2c5-15b19d7116d5","type":"LinearAxis"},{"id":"97ec3a3c-01c7-4912-8e62-98284d63dae3","type":"Grid"},{"id":"9be237e4-696b-402b-a816-cca962343a6f","type":"BoxAnnotation"},{"id":"3bdbbc32-717b-4a00-bf0d-834f224c5daf","type":"Legend"},{"id":"9512ea20-4566-45a9-aa77-6b3a170237c9","type":"GlyphRenderer"},{"id":"0753661f-95aa-4388-b6bc-2de50ae27dd9","type":"GlyphRenderer"},{"id":"6d8ec9f4-f4f2-4590-94e2-4170ef2cd166","type":"GlyphRenderer"},{"id":"4088aadd-2cbd-466f-a523-a9c87a1566b0","type":"LabelSet"}],"title":{"id":"9135a977-2a41-4547-92f8-3a2959da633a","type":"Title"},"toolbar":{"id":"6fd258fe-34c2-43e3-9581-aab3946e21d0","type":"Toolbar"},"x_range":{"id":"d3d39b15-11c0-4760-9b8a-5eee27de6a51","type":"DataRange1d"},"x_scale":{"id":"7090e1e7-53a1-4532-b87c-5e8e94868d5c","type":"LinearScale"},"y_range":{"id":"13a72d3f-fabe-419c-978c-f00870aef8b0","type":"DataRange1d"},"y_scale":{"id":"cd033eeb-4342-4644-8a81-6e4a57836b62","type":"LinearScale"}},"id":"8f41f3fb-4ac9-438c-9fc4-526047dfe786","subtype":"Figure","type":"Plot"},{"attributes":{},"id":"7090e1e7-53a1-4532-b87c-5e8e94868d5c","type":"LinearScale"},{"attributes":{},"id":"55962ce4-33ed-421b-85e9-c1b7b72fdd2a","type":"BasicTicker"},{"attributes":{},"id":"96fd0775-2d63-4d36-9dcc-f8c3f0ca0ca0","type":"ResetTool"},{"attributes":{},"id":"702fd5c2-b7c4-46f4-ac99-a654b7a61a1e","type":"BasicTicker"},{"attributes":{"overlay":{"id":"9be237e4-696b-402b-a816-cca962343a6f","type":"BoxAnnotation"}},"id":"cb5762ed-e9a4-4774-94d2-d66eca83117c","type":"BoxZoomTool"},{"attributes":{"source":{"id":"14a6fc3b-7f92-488c-81d4-fe802b1b1b15","type":"ColumnDataSource"}},"id":"5cb07060-2988-4b92-b7a5-a3595b59e9e7","type":"CDSView"},{"attributes":{"axis_label":"Partial Pressure CO2 (Pa)","formatter":{"id":"ea1bf303-86c0-4fe6-bce1-70f28f526db1","type":"BasicTickFormatter"},"plot":{"id":"8f41f3fb-4ac9-438c-9fc4-526047dfe786","subtype":"Figure","type":"Plot"},"ticker":{"id":"702fd5c2-b7c4-46f4-ac99-a654b7a61a1e","type":"BasicTicker"}},"id":"40860ef2-c2cc-40ed-a2c5-15b19d7116d5","type":"LinearAxis"},{"attributes":{"label":{"value":"Bulk vapor"},"renderers":[{"id":"9512ea20-4566-45a9-aa77-6b3a170237c9","type":"GlyphRenderer"}]},"id":"27dc17d5-2d1c-4db5-a830-edf5679012e2","type":"LegendItem"},{"attributes":{"axis_label":"Axial distance from top (m)","formatter":{"id":"bcc2034d-131e-4153-b7b9-d98f10a07037","type":"BasicTickFormatter"},"plot":{"id":"8f41f3fb-4ac9-438c-9fc4-526047dfe786","subtype":"Figure","type":"Plot"},"ticker":{"id":"55962ce4-33ed-421b-85e9-c1b7b72fdd2a","type":"BasicTicker"}},"id":"18ba8f4c-a5ea-477a-931e-d58496934361","type":"LinearAxis"},{"attributes":{"callback":null,"column_names":["y","x","labels"],"data":{"labels":["Reboiler vapor"],"x":[12.0956],"y":[10827.3815090278]}},"id":"8346089e-d296-4764-91ec-ff4ef05a917b","type":"ColumnDataSource"},{"attributes":{"data_source":{"id":"c4414caf-b280-422b-8483-1807c1d68dfe","type":"ColumnDataSource"},"glyph":{"id":"7bcce2c0-0680-4927-b068-e68d78fd7202","type":"Line"},"hover_glyph":null,"muted_glyph":null,"nonselection_glyph":{"id":"f02a2f3f-3c10-4e60-b23e-3cddc2528ddf","type":"Line"},"selection_glyph":null,"view":{"id":"795116a0-8abf-46eb-99bf-a882e1d1b744","type":"CDSView"}},"id":"0753661f-95aa-4388-b6bc-2de50ae27dd9","type":"GlyphRenderer"},{"attributes":{"callback":null,"column_names":["x","y"],"data":{"x":{"__ndarray__":"AAAAAAAAAAB8/stb1DPBP3z+y1vUM9E/u/2xib7N2T98/stb1DPhPyWZhdDigOU/xZh459fN6T9jmGv+zBruPwFMrwrhM/E/0csollta8z+gS6Ih1oD1P3DLG61Qp/c/QEuVOMvN+T8Oyw7ERfT7P2OYa/7MGv4/GozyxKMgAEABTK8K4TMBQOkLbFAeRwJA0csolltaA0C4i+XbmG0EQKBLoiHWgAVAiAtfZxOUBkBwyxutUKcHQBoySkqUughAAvIGkNHNCUDqscPVDuEKQNJxgBtM9AtAuTE9YYkHDUCi8fmmxhoOQImxtuwDLg9AuLg5maAgEECsGBg8P6oQQAFMrwrhMxFA9auNrX+9EUDpC2xQHkcSQN1rSvO80BJA0csolltaE0DFKwc5+uMTQLiL5duYbRRArevDfjf3FECgS6Ih1oAVQPZ+OfB3ChZA6d4XkxaUFkDdPvY1tR0XQNGe1NhTpxdAxf6ye/IwGEC5XpEekboYQKy+b8EvRBlAoB5OZM7NGUCVfiwHbVcaQOqxw9UO4RpA3hGieK1qG0DScYAbTPQbQMXRXr7qfRxAuTE9YYkHHUCtkRsEKJEdQKLx+abGGh5AlVHYSWWkHkCJsbbsAy4fQN7kTbultx9AaSIWL6IgIEBjUoWAcWUgQF2C9NFAqiBAVrJjIxDvIEBR4tJ03zMhQEsSQsaueCFARUKxF369IUA+ciBpTQIiQOkLbFAeRyJA4zvboe2LIkDda0rzvNAiQNebuUSMFSNA0csolltaI0DL+5fnKp8jQMUrBzn64yNAv1t2iskoJEC4i+XbmG0kQGMlMcNpsiRAXVWgFDn3JEBXhQ9mCDwlQFG1frfXgCVAS+XtCKfFJUBFFV1adgomQD9FzKtFTyZAOXU7/RSUJkAypapO5NgmQN0+9jW1HSdA2G5lh4RiJ0DRntTYU6cnQMvOQyoj7CdAxf6ye/IwKEA=","dtype":"float64","shape":[91]},"y":{"__ndarray__":"tPRumMz54UB+g/DBIQnmQD04nW2/oehAfJJ3Qd8v6kAsqlQslxXrQBJmdERdl+tA8SDmX8zf60AS7nfW9AfsQPMoAbogHuxAi8BboVQq7EA3H70LBjHsQL5c6KOsNOxAPqdR8aQ27ECvdCDhrjfsQNJoFvIzOOxAMnzX92047ED6o6CAfDjsQNjuppdwOOxAJJpzQFM47ED8z9oFKTjsQH64Wu3zN+xAf+qmh7Q37EBrjQyEajfsQD8WO/oUN+xAq6RopLI27ED+cULTQTbsQCB7e4nANexAIpBKcyw17EBRu33fgjTsQHK8/bTAM+xAkd9hZeIy7ECMj8fc4zHsQD+FM2DAMOxABS7btXIv7EB+thWn9C3sQD1/oi0/LOxAWXPHPEoq7EBh0UqbDCjsQF9+0Ld7JexAWR7ldosi7EAAiOr5LR/sQOLdeS1TG+xA2IAfTekW7ECzhbZC2xHsQCABIT0RDOxAM+8I/28F7EC7MEBp2P3rQBIyt/gm9etA0My+ODPr60As3jsrz9/rQDcBWBDG0utAY5s9Tt7D60BtbZS+1LLrQHfGfyRfn+tAwDz13iqJ60A730Wr3G/rQGFgp6sQU+tAzpvUyloy60BlIECXRw3rQMntsAdc4+pATQBvgR606kCZ8FlyDn/qQO98N+2vQ+pAh2w9Vo0B6kB/8qjfO7jpQDRdETFgZ+lAOogy6bIO6UDVOZ2fBK7oQHAabrc9RehA677FTG7U50DrtEyxuFvnQGcbzLxg2+ZAmzc1/MRT5kCrGyJgW8XlQIMeBvWsMOVA6W3P/FCW5EDhEJnC5/bjQMUNyBsRU+NAOjBqt3mr4kBjGd8dugDiQNvbfJRnU+FA/F1VEguk4ECfgoiePebfQC1AtmkYgt5AHCdd31Mc3UDY2LMOc7XbQH1cxKK9TdpA67MWi1Dl2ECfyd3zuXvXQEkx3Ee0D9ZAaeGG3hae1EA=","dtype":"float64","shape":[91]}}},"id":"c4414caf-b280-422b-8483-1807c1d68dfe","type":"ColumnDataSource"}],"root_ids":["8f41f3fb-4ac9-438c-9fc4-526047dfe786"]},"title":"Bokeh Application","version":"0.12.9"}};
                      var render_items = [{"docid":"80139986-1908-4271-886f-f4e38b644bfb","elementid":"e3811af1-a8f9-4042-b961-7400b192ed89","modelid":"8f41f3fb-4ac9-438c-9fc4-526047dfe786"}];

                      root.Bokeh.embed.embed_items(docs_json, render_items);
                    }

                    if (root.Bokeh !== undefined) {
                      embed_document(root);
                    } else {
                      var attempts = 0;
                      var timer = setInterval(function(root) {
                        if (root.Bokeh !== undefined) {
                          embed_document(root);
                          clearInterval(timer);
                        }
                        attempts++;
                        if (attempts > 100) {
                          console.log("Bokeh: ERROR: Unable to embed document because BokehJS library is missing")
                          clearInterval(timer);
                        }
                      }, 10, root)
                    }
                  })(window);
                });
              };
              if (document.readyState != "loading") fn();
              else document.addEventListener("DOMContentLoaded", fn);
            })();

            </script>
        </body>
    </html>