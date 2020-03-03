" Location ~/.config/nvim/init.vim

set runtimepath^=~/.vim runtimepath+=~/.vim/after
let &packpath = &runtimepath

"----------------------------
" START PLUGINS
" install vim-plug
if empty(glob('~/.vim/autoload/plug.vim'))
    silent !curl -fLo ~/.vim/autoload/plug.vim --create-dirs
      \ https://raw.githubusercontent.com/junegunn/vim-plug/master/plug.vim
    autocmd VimEnter * PlugInstall --sync | source ~/.vimrc
endif

if empty(glob('~/.local/share/nvim/site/autoload'))
    silent !curl -fLo ~/.local/share/nvim/site/autoload --create-dirs
      \ https://raw.githubusercontent.com/junegunn/vim-plug/master/plug.vim
    autocmd VimEnter * PlugInstall --sync | source ~/.vimrc
endif



call plug#begin('~/.vim/plugged')

" R-related 
Plug 'jalvesaq/Nvim-R'
Plug 'roxma/nvim-yarp'
Plug 'ncm2/ncm2'
Plug 'vim-pandoc/vim-pandoc-syntax'

"This the original one
"Plug 'gaalcaras/ncm-R'
"Corrects for the anoying massage in case nvim hasn't started
Plug 'ShawnChen1996/ncm-R'

Plug 'sirver/UltiSnips'
Plug 'ncm2/ncm2-ultisnips'
Plug 'pablo-gar/vim-snippets'
Plug 'lervag/vimtex'

" Fast python completion (use ncm2 if you want type info or snippet support)
Plug 'HansPinckaers/ncm2-jedi'

" Nice status bar
Plug 'vim-airline/vim-airline'
Plug 'vim-airline/vim-airline-themes'

" Words in buffer completion
Plug 'ncm2/ncm2-bufword'
" Filepath completion
Plug 'ncm2/ncm2-path'

" Syntax checking
"Plug 'dense-analysis/ale'

" Map shift+j for ultsnips (snippets)
let g:UltiSnipsExpandTrigger = "<C-j>"
let g:UltiSnipsJumpForwardTrigger="<C-j>"
let g:UltiSnipsJumpBackwardTrigger="<C-k>"
let g:UltiSnipsSnippetDirectories=["UltiSnips", $HOME."/.vim/plugged/my_snippets/"]

" Initialize plugin system
call plug#end()



" END PLUGINS
"----------------------------
"
let g:python3_host_prog="/home/users/paedugar/.conda/envs/pybase/bin/python3"
colorscheme peachpuff
    

"-----------------------------
" ncm2 CONFIG
autocmd BufEnter * call ncm2#enable_for_buffer()
set completeopt=noinsert,menuone,noselect
" make it fast
set shortmess+=c
let ncm2#popup_delay = 5
let ncm2#complete_length = [[1, 1]]
set pumheight=5
"Map tab for autocompletion
inoremap <expr> <tab> pumvisible() ? "\<C-n>" : "\<Tab>"
inoremap <expr> <S-tab> pumvisible() ? "\<C-p>" : "\<S-Tab>"
inoremap <silent> <expr> <CR> (pumvisible() && empty(v:completed_item)) ?  "\<c-y>\<cr>" : "\<CR>"

"-----------------------------
" Re-mapping key bindings
"remap esc to jk in terminal
tnoremap jk <C-\><C-n>

" general
inoremap (( ()<ESC>i
inoremap [[ []<ESC>i
inoremap '' ''<ESC>i
inoremap "" ""<ESC>i
inoremap jk <ESC>
" insert date
nnoremap <F2> a<C-R>=strftime("%c")<CR><Esc>
inoremap <F2> <C-R>=strftime("%c")<CR><Esc>
let mapleader = "\<Space>"
set number
set hlsearch

"-----------------------------
" Syntax 

"Markdown and Rmarkdown
augroup pandoc_syntax
    au! BufNewFile,BufFilePre,BufRead *.md set filetype=markdown.pandoc
augroup END
 
" Snakemake highlighting syntax
autocmd BufNewFile,BufRead Snakefile set syntax=snakemake
autocmd BufNewFile,BufRead *.smk set syntax=snakemake

"-----------------------------
" General

"syntax on
set encoding=utf-8
set undodir=~/.vim/undodir
set undofile  " save undos
set undolevels=10000  " maximum number of changes that can be undone
set undoreload=100000  " maximum number lines to save for undo on a buffer reload

" enable filetype detection:
filetype plugin indent on

"split navigations
"nnoremap <C-J> <C-W><C-J>
"nnoremap <C-K> <C-W><C-K>
"nnoremap <C-L> <C-W><C-L>
"nnoremap <C-H> <C-W><C-H>

" Allows faster switching between windows with "Ctrl+w"! By Questor
nnoremap <silent> <C-w> <C-w><C-w>
inoremap <silent> <C-w> <Esc><C-w><C-w>
vnoremap <silent> <C-w> <Esc><C-w><C-w>

" opened where left off 
if has("autocmd")
  au BufReadPost * if line("'\"") > 1 && line("'\"") <= line("$") | exe "normal! g'\"" | endif
endif

" Insert new line while keeping identation
inoremap <CR> <CR>x<BS>
nnoremap o ox<BS>
nnoremap O Ox<BS>

" <Ctrl-l> redraws the screen and removes any search highlighting.
nnoremap <silent> <C-l> :nohl<CR><C-l>

" Tab == n spaces
set expandtab
set tabstop=4
set shiftwidth=4

" Make backsapce key work
set backspace=2

" force tabs in make files
autocmd FileType make set noexpandtab
autocmd FileType makefile set noexpandtab



"----------------------------
" Nvim-R config
 
" remapping the basic :: send line
nmap <Space> <Plug>RDSendLine

" remapping selection :: send multiple lines
vmap <Space> <Plug>RDSendSelection

" remapping selection :: send multiple lines + echo lines
vmap <Space>e <Plug>RESendSelection

" remapping sending a markdonw chunk
nmap <C-Space> <Plug>RSendChunk

" To go into normal mode while in R terminal do <C-\><C-N>
let R_esc_term = 0

" Open R terminal in tmux buffer (0). Currently does not work, you have to
" install a terminal emulator
let R_in_buffer = 1

" Clears line before executing a new one
let R_clear_line = 1

" Disables automatic <- when typing _
let R_assign = 0

" Do horizontal split
let R_rconsole_width = 0
let R_rconsole_height = 20

" END Nvim-R config
"----------------------------

"---expandt-----
" ale config (syntax correction)

" Enable nice status through airline
"let g:airline#extensions#ale#enabled = 1

"Move between errors mapping
nmap <silent> <C-k> <Plug>(ale_previous_wrap)
nmap <silent> <C-j> <Plug>(ale_next_wrap)

" Check Python files with flake8 and pylint.
let b:ale_linters = ['autopep8', 'isort', 'yapf']
" Fix Python files with autopep8 and yapf.
let b:ale_fixers = ['autopep8', 'isort', 'yapf']

" END Ale config (syntax correction)
"----------------------------


"----------------------------
" vim-airline config (status bar)

let g:airline_theme='papercolor'

" END vim-airline config (status bar)
"---expand
