[Back to start](./README.md)


Provided Tools
--------------

You'll find Go packages in subfolders of this repo:

* `src/genomics/common/` Data models, and common routines across all tools.
* `src/genomics/bed/` [Bed docs](./bed.md)
* `src/genomics/align/` [Align docs](./align.md)
* `src/genomics/reduce/` [Reduce docs](./reduce.md)
* `src/genomics/vc/` [VC docs](./vc.md)

Building and Running individual tools
---------------------------------------

1. heads up: code's written in go.  The [Helpful Links](#helpful-links) section covers some essentials.

1. Get a recent version of golang (1.9.x):

    ```bash
    sudo add-apt-repository ppa:gophers/archive
    sudo apt update
    sudo apt-get install golang-1.9-go

    # that above installs it in /usr/lib/go-1.9/bin
    # you can update ~/.bashrc, and update your current shells
    # on the spot:
    export PATH=/usr/lib/go-1.9/bin:"$PATH"
    ```

    Make sure you have the right version:

    ```bash
    $ go version
    go version go1.9.2 linux/amd6
    ```

1. For newcomers to golang, make sure your `$GOPATH` is setup adequately. (e.g. there should exist a directory DIR in your $GOPATH, such that $DIR/src/genomics is the root of this repo). 

   If you're new to go, it's likely to surprise you. You may read up on gopath: `go help gopath`

   If you're used to having all the repos you contribute to under a
   common folder, regardless of language (e.g. `~/src/` containing
   Java, C, JavaScript, Python, etc. programs), you will appreciate
   configuring more than one directory in your `$GOPATH`
   (e.g. `GOPATH=DIR1:DIR2`).  The first DIR in `$GOPATH`, by default,
   is the one where dependencies downloaded with `go get` will end up
   -- this one directory (DIR/src/) will get invariably "polluted"
   with tons of github repos and various cloud-downloaded packages.

   When building, running, or importing packages, every `DIR/src/`
   will be consulted.  If some of these packages are installed,
   they'll end up in the respective `DIR/bin/`, so, make sure you
   update your `$PATH` in a way that reflects all DIRs in
   your`$GOPATH` directories.

   For instance, if you have `GOPATH="DIR1:DIR2"`, you will want:

   ```bash
   export PATH="DIR1/bin:DIR2/bin:$PATH"
   ```

1. compile one of the packages

    ```bash
    cd align  # or any other package
    go get    # will download deps
    go build  # should leave you with a binary in the folder.
    ```

1. Run a tool:

    ```bash
    ./align --help
    ```

Helpful links
---------

* The go way, written by a programmer [https://github.com/alco/gostart]
* The go way, written by golang [https://golang.org/doc/code.html]


Setting up your editor
---------------------

There is a generic list of IDEs and Plugins for go here: https://github.com/golang/go/wiki/IDEsAndTextEditorPlugins

### Emacs

* This is a pretty good summary : [http://andrewjamesjohnson.com/configuring-emacs-for-go-development/]


  The code below is a summary of the above, minus the bugs and typos.

  ```elisp

  ;; this assumes you have the MELPA package module setup
  (package-initialize)
  ...

  ;; Youll have to run the following outside emacs
  ;;   go get -u -v github.com/rogpeppe/godef
  ;;   go get -u -v github.com/nsf/gocode
  ;;   go get -u -v golang.org/x/tools/cmd/godoc

  ;; You'll have to run the following in emacs
  ;;   M-x package-install go-mode (installs go-mode and godef)
  ;;   M-x package-install go-eldoc
  ;;   M-x package-install auto-complete
  ;;   M-x package-install go-autocomplete


  (defun go-mode-setup ()
    (setq compile-command "go build -v && go test -v && go vet")
    (define-key (current-local-map) "\C-c\C-c" 'compile) ;; Ctrl-c Ctrl-c to compile
    (go-eldoc-setup)                                     ;; show argument of function at point
    (add-hook 'before-save-hook 'gofmt-before-save)      ;; run gofmt on save
    (local-set-key (kbd "M-.") 'godef-jump))             ;; M-. to jump to symbol. M-* to return.

  (add-hook 'go-mode-hook 'go-mode-setup)

  (require 'auto-complete-config)
  (require 'go-autocomplete)
  (ac-config-default)
  ```
  For live error checking, you should consider installing flycheck which has good go support:

  ```elisp
  (use-package flycheck
    :ensure t)
  (add-hook 'after-init-hook #'global-flycheck-mode)
  ```
  
  Once you're in a buffer in go-mode, you may want to check that the
  setup is correct with `M-x flycheck-verify-setup`, which will print
  the paths used for each checker found on your system. You may want to
  install the following checkers:
  
  * Go Lint: `go get -u github.com/golang/lint/golint`
  * Go errcheck (ensure error checking happens): `go get -u github.com/kisielk/errcheck`
  * Go Unconvert (identify unnecessary type conversions): `go get -u github.com/mdempsky/unconvert`
  * Go megacheck (staticheck + gosimple + unused): `go get -u honnef.co/go/tools/cmd/megacheck`
