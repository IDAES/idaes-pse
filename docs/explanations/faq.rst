Frequently Asked Questions
==========================

How do I...
-----------

... Run some examples?
    See :doc:`../tutorials/tutorials_examples`.

... Get more help?
    Use the website to `register <https://idaes.org/register/>`_ for the IDAES support mailing list.
    Then you can send questions to idaes-support@idaes.org. For more specific technical questions, we recommend
    our new `IDAES discussions board on Github <https://github.com/IDAES/idaes-pse/discussions>`_.

Troubleshooting
---------------

Missing win32api DLL
    For Python 3.8 and maybe others, you can get an error when running Jupyter on Windows 10 about
    missing the win32api DLL. There is a relatively easy fix::

        pip uninstall pywin32
        pip install pywin32==225
