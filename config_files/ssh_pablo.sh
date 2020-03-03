# Location ~/.ssh/config
Host sher2
    HostName sh-ln01.stanford.edu
    User paedugar
    ForwardX11 yes
    ServerAliveInterval 60
    ForwardAgent yes
    GSSAPIAuthentication yes
    GSSAPIDelegateCredentials yes
    ControlMaster auto
    ControlPath ~/.ssh/%r@%h:%p
    ControlPersist yes

Host sher_dtn sherlock2_dtn dtn.sherlock.stanford.edu
    HostName dtn.sherlock.stanford.edu
    User paedugar
    ForwardX11 no
    ServerAliveInterval 60
    ForwardAgent yes
    GSSAPIAuthentication yes
    GSSAPIDelegateCredentials yes
    ControlMaster auto
    ControlPath ~/.ssh/%r@%h:%p
    ControlPersist yes

Host sher
    HostName login.sherlock.stanford.edu
    User paedugar
    ForwardX11 yes
    ServerAliveInterval 60
    ForwardAgent yes
    GSSAPIAuthentication yes
    GSSAPIDelegateCredentials yes
    ControlMaster auto
    ControlPath ~/.ssh/%r@%h:%p
    ControlPersist yes

Host *
  AddKeysToAgent yes
  UseKeychain yes
  IdentityFile ~/.ssh/id_rsa
