# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/opt/software/anaconda/python-3.9.7/2021.11/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/opt/software/anaconda/python-3.9.7/2021.11/etc/profile.d/conda.sh" ]; then
        . "/opt/software/anaconda/python-3.9.7/2021.11/etc/profile.d/conda.sh"
    else
        export PATH="/opt/software/anaconda/python-3.9.7/2021.11/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

