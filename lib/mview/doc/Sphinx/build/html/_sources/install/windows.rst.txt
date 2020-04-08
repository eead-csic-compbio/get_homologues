Installation: Windows
=====================

There are two ways to install MView.

Either method can be used by an ordinary user installing into their own
account, or by a system administrator installing onto a computer with multiple
users. It is assumed that Perl is already installed and on your ``PATH``.

* :ref:`Installer script<ref_windows_installer>`
* :ref:`Manual install<ref_windows_manual>`
* :ref:`How to set PATH<ref_windows_path>`


.. _ref_windows_installer:

Installer script
^^^^^^^^^^^^^^^^

The installer program should work on all systems, but is new and relatively
experimental.

You unpack the archive into a destination folder and run the installer from
there, following the instructions. You may have to edit ``PATH`` afterwards.

Explanation: the installer puts a small mview driver program into a folder on
``PATH`` so that it can be run easily by the user. The driver knows the
location of the unpacked MView folder and starts the real MView program.

1. Save the archive to somewhere under your home folder then uncompress and
   extract it (using an archiver like 7-Zip, as here)::

        7z x mview-VERSION.zip

   This creates a sub-folder ``mview-VERSION`` containing all the files.
   
2. Change to this folder.

3. Run the command::

        perl install.pl
        
   and follow the instructions. You will be offered various places to install
   the driver script.
   
   If you know in advance the name of the folder you want to use for the
   driver script, you can supply it on the command line::

        perl install.pl drive:\folder\on\my\path

4. If the installer couldn't find a sensible place to install the driver, it
   chooses ``C:\bin`` and you will have to add that to your ``PATH``, then
   start a new command prompt.


.. _ref_windows_manual:

Manual install
^^^^^^^^^^^^^^

1. Save the archive to your software area, for example, ``C:\Program Files``,
   then uncompress and extract it (using an archiver like 7-Zip, as here)::

        7z x mview-VERSION.zip

   This creates a sub-folder ``mview-VERSION`` containing all the files.

2. Change to this folder.

3. Edit the file ``bin\mview``.

  * Find the line::

        $MVIEW_HOME = "/path/to/mview/unpacked/folder";

    and change the path, in our example, to::

        $MVIEW_HOME = "C:\Program Files\mview-VERSION";

  * Save the file.

4. Finally, make sure that the ``bin`` folder containing the ``mview`` script
   (that you just edited) is on the user ``PATH``, then start a new command
   prompt.

   In our example, you would add ``C:\Program Files\mview-VERSION\bin`` to the
   existing value of ``PATH``, or replace any older MView path.


.. _ref_windows_path:

How to set PATH
^^^^^^^^^^^^^^^

The ``PATH`` environment variable is a list of ``;`` (semi-colon) separated
folders containing programs. When you type the name of a program at the
command prompt, the system searches these folders, in order, until it finds
the program and runs it (or complains if the program can't be found).

Assume you are adding ``C:\bin`` as the directory containing the newly
installed mview script. On all systems the ``PATH`` environment variable would
be extended by adding ``C:\bin`` to the existing ``PATH`` value using
semi-colon delimiters as needed. You can prepend the new path (it will be
searched first for commands), insert it somewhere in the middle, or append it
at the back (it will be searched last).

How to change ``PATH`` on different editions of Windows (gleaned from the
Internet):

* Windows 10 and Windows 8

  - Search » System (Control Panel) » Advanced system settings »
    Environment Variables

* Windows 7

  - Computer » Properties » Advanced system settings » Environment
    Variables

* Windows Vista

  - My Computer » Properties » Advanced system settings » Environment
    Variables

* Windows XP

  - Start » Control Panel » System » Advanced » Environment Variables

On all systems, once you've updated the ``PATH`` variable, open a new command
prompt, then the ``mview`` command should be recognised, so that running::

  mview -help

prints the help message for the new version.

Note: if you already have an older mview installed on the ``PATH`` and append
the new location at the back of ``PATH``, the older program will still be
found first whenever you try to run mview, so be aware of that; you would need
to delete the old version, or rearrange the ``PATH`` order.
