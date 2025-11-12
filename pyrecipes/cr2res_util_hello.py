from typing import Any, Dict

# Import the required PyCPL modules
import cpl.ui


# Define our "Hello, world!" recipe as a class which inherits from
# the PyCPL class cpl.ui.PyRecipe
class HelloWorld(cpl.ui.PyRecipe):
    # The information about the recipe needs to be set. The base class
    # cpl.ui.PyRecipe provides the class variables to be set.
    # The recipe name must be unique, because it is this name which is
    # used to identify a particular recipe among all installed recipes.
    # The name of the python source file where this class is defined
    # is not at all used in this context.
    _name = "helloworld"
    _version = "1.0"
    _author = "U.N. Owen"
    _email = "unowen@somewhere.net"
    _copyright = "GPL-3.0-or-later"
    _synopsis = "PyCPL version of 'hello, world!'"
    _description = (
        "This is the PyCPL version of the well known hello, world program.\n"
        + "It says hello to each input file in the input set-of-frames."
    )

    # Our recipe class also needs to provide the run() method with the
    # correct arguments and return values.
    #
    # As inputs the run method must accept a cpl.ui.FrameSet, and a dictionary
    # that contains the parameters of the recipe.
    # In this example the rcipe does not have any recipe parameters, but the
    # function still has to accept the dictionary as second argument.
    #
    # When the recipe is done, it has to return the produced product files in
    # another cpl.ui.FrameSet object. The current example will not create any
    # output data, so that it will have to return an empty cpl.ui.FrameSet.
    def run(
        self, frameset: cpl.ui.FrameSet, settings: Dict[str, Any]
    ) -> cpl.ui.FrameSet:
        for frame in frameset:
            print(f"Hello, {frame.file}!")
        return cpl.ui.FrameSet()
