#!/bin/bash
# clang-format does not indent #pragmas (or other #preprocessor directives) with the code.
# This small script indents pragmas (BEWARE: very long pragma args wrap if clang-format 
# reflows comments!)

sed --expression='s/#pragma/\/\/ #pragma/' | clang-format | sed --expression='s/\/\/ #pragma/#pragma/'

# This works by comments:  #pragma -> // #pragma -> #pragma
# which works because comments are properly aligned.
#
# For vscode, enable with Custom Local Formatter -extension, with json settings
#   "customLocalFormatters.formatters": [
#        {        
#            "command": "/path/to/clang-format-with-pragmas.sh",
#            "languages": ["cpp", "c"]
#        }
#    ],
#    "[cpp]": {
#        "editor.defaultFormatter": "jkillian.custom-local-formatters"
#    }
