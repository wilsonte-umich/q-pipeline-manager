#!/usr/bin/perl
use strict;
use warnings;
$ENV{PATH} = "/garage/wilsonte_lab/bin/q/q-pipeline-manager-develop/utilities/bedutil-1.0.0:/garage/wilsonte_lab/bin/q/q-pipeline-manager-develop/utilities/slurp-1.0.0:/garage/wilsonte_lab/bin/q/q-pipeline-manager-develop/utilities/stat-1.0.0:$ENV{PATH}";
our $version = '1.0.5';
our $perlPath = '/usr/bin/perl';
our $bashPath = '/usr/bin/bash';
our $timePath = '/usr/bin/time';
our $timeVersion = '1.8';
our $memoryCorrection = '1';
our $memoryMessage = '';
our $qDir = '/garage/wilsonte_lab/bin/q/q-pipeline-manager-develop';
our $libDir = '/garage/wilsonte_lab/bin/q/q-pipeline-manager-develop/lib';
our $modulesDir = '/garage/wilsonte_lab/bin/q/q-pipeline-manager-develop/modules';
our $qType = 'SGE';
our $schedulerDir = '/opt/sge/bin/lx-amd64';
$ENV{Q_Q_TYPE} = 'SGE';
$ENV{Q_UTIL_DIR} = '/garage/wilsonte_lab/bin/q/q-pipeline-manager-develop/utilities';
$ENV{Q_MOD_DIR} = '/garage/wilsonte_lab/bin/q/q-pipeline-manager-develop/modules';
our %shellCommands = map { $_ => 1 } qw(cpan sqldiff netctl-auto sha512sum kernel-install whatis p11tool groupmems nf-log postcat rst2html5 smbd filtersubs.sh fusermount lzmainfo lsmem pcprofiledump ngettext systemd-umount xargs hexdump snap_scheduler.py unicode_start numademo pvresize x86_64-pc-linux-gnu-gcc-nm pyrsa-decrypt telnet gpgrt-config pkgdelta gnutls-serv config kmercountmulti.sh pod2text vde_pcapplug net makeconv nl-monitor bind9-config callgenes.sh loginctl rpcdebug pinentry-tty man mdb_load convert-mans nf-queue sndfile-deinterleave deactivate_matlab.sh STAR x86_64-pc-linux-gnu-c++ xfs_quota grub-mkimage samba_dnsupdate guild autrace fc-validate iptables-nft-save qrls getreads.sh systemd-cgls grodvi cifsiostat watchgnupg md5sum xtotroff perlthanks bowtie2-align-s-debug readlength.sh reducesilva.sh shrinkaccession.sh ftp pdfmom giffilter locale processfrag.sh pg_dump get_module telnetd yapp lvresize nl-neigh-delete colrm zcmp filterbytaxa.sh nl-link-name2ifindex pluginviewer ldns-nsec3-hash arpaname h2xs nl-class-add gss-client eqn2graph db_archive mergesketch.sh nl-link-stats gettext.sh ssltap muxbyname.sh filecap ldapmodrdn pathchk nmap testfilesystem.sh create-cracklib-dict ptargrep logname filtersam.sh rpcbind streamsam.sh dbus-run-session era_restore sndfile-convert ausyscall daisync-plot sfdisk lscpu grotty mailq rst2man gfortran pl2pm klist cairo-trace bamToBed tmux gdbus mysqladmin xmbind chsh testformat.sh tiffcrop tar netctl pkaction pg_restore update-ca-trust shlibsign pigz gc-ctrl.py display-coords rst2s5 glib-mkenums insmod id killall ping mdadm showwal dirfingerprint ksu q pk12util systemd-socket-activate arping imgcmp nsupdate usbredirserver shuffle.sh qmqp-sink link db_printlog tac dig saslpasswd2 gcov runlevel bbsketch.sh neqn gdbus-codegen jpgicc c++ pacman unzip.sh genrb bowtie2-align-l-debug gcc-ranlib iptables-nft true lvreduce ldapmodify lvmpolld slapauth od pkgconf modprobe nl-classid-lookup lvrename osd_login ssh-keygen journalctl grub-install sensors-detect phylip2fasta.sh locktest postsuper x86_64 vde_autolink lz4c pacman-conf showconsolefont cksum iconvconfig countbarcodes.sh bbversion.sh xfs_estimate e4crypt clear gpg2 complementBed lvmdiskscan pdfroff ldns-verify-zone configure.pl rmdir e2scrub_all ptar pkttyagent glib-genmarshal grub-bios-setup replaceheaders.sh comparegff.sh bowtie2-inspect-l-debug thin_rmap python3.7 setup.py gpasswd kdestroy pvck gawk-4.2.1 flock svtools mv closestBed pskctool thin_trim bbmapskimmer.sh genomeCoverageBed screen pg_dumpall disable-paste fancontrol hostname set_stat kbxutil dmidecode nl-cls-delete idn2 c_rehash sudoreplay loadunimap bsdtar nettle-hash hwloc-patch qdel ifstat sg factor mvxattr pslog tty uncompress rndc-confgen setpci plipconfig qmon gd2copypal mail sexp-conv fc-pattern slencheck migratepages colcrt nl-link-ifindex2name sha384sum fastaFromBed a_sample_mt.sh diff bowtie2-inspect-l nl-route-delete sge_execd gd2togif linux64 pwdx libgcrypt-config hostnamectl dltest glilypond ntlm_auth grub-file bridge gencat grep lzfgrep kadmind vgmerge psfgettable ntpdate arptables-nft-save scp msvtools zegrep json_pp COMPILE nss-config ldns-testns vdeterm faillog setpriv gcc-nm smbtar ul more sotruss repo-remove healthd zipnote sftp mount.cifs pg_config mkdir fax2ps analyzeaccession.sh vgchange postkick shred.sh loglog.sh mouse-test grub-mkrescue json-glib-format vgs bbsplitpairs.sh ntpdc mii-tool dropuser rpc.statd gtk-launch splitbytaxa.sh zipinfo python3.7-config gf_attach ip6tables-translate script choom elfedit vgcreate qrsh rqcfilter2.sh networkctl info xfs_logprint zfgrep devlink iptunnel sge_qmaster zcat python2-config newuidmap pod2man rst2xml arptables-nft-restore gdk-pixbuf-query-loaders winbindd dconf genl fc-cache pyvenv-3.7 getcifsacl chvt postconf swapoff pairToBed dumpsexp dnssec-keygen diskbench.sh zstdgrep setarch coverageBed kadmin fribidi smooth xzcat ldbsearch pyrsa-encrypt fuser nl-addr-add avahi-browse-domains analyzegenes.sh nl-class-list getpcaps texindex gencfu ldns-zsplit pr catman crlutil uptime sendmail msgcmp dbilogstrip rst2pseudoxml gtk-builder-tool smbtree js52 paplay mysql bowtie2-build-s sha224sum celtenc051 dumpe2fs makechimeras.sh basename nf-ct-add calctruequality.sh dd gdk-pixbuf-pixdata vgexport ldd mergepgm.sh unixterm filtervcf.sh systemd-notify gluster-eventsapi nl-addr-delete kmerlimit.sh jemalloc-config msgen lz4cat umount.nfs4 lspci pager rview iconv json-glib-validate hb-ot-shape-closure fdisk countsharedlines.sh zenmap sndfile-metadata-get pic2graph avahi-resolve ssh taxsize.sh xzfgrep pwck pyrsa-keygen rlogin-cwd gettextize vgcfgbackup vpddecode qquota addadapters.sh qsh fc-conflist view grub-mount mke2fs cryptdir debugreiserfs psfstriptable c++filt rst2odt_prepstyles pngfix ldns-walk gssproxy zipcloak summarizecontam.sh resize_reiserfs msgmerge mdig bbstats.sh pvdisplay dwp thin_delta nfs-ls reindexdb rfkill webpmux gluster-mountbroker env_parallel nspr-config mbni gdk-pixbuf-thumbnailer cp tagBam tail env_parallel.tcsh samba e4defrag ln makepkg-template sim_server kadmin.local guile-tools hwclock grog curl-config broadwayd pre-grohtml avahi-set-host-name autopoint dnssec-verify bowtie2-inspect idle3 update-smart-drivedb aulastlog refer xfs_rtcp git-shell gio-launch-desktop tiff2pdf mkinstalldirs lvdisplay csslint-0.6 grub-menulst2cfg envsubst grub-macbless chcon filterassemblysummary.sh ip6tables tiffdump timed-run gdbm_load busy-nodes lsusb.py mkswap setsid nf-exp-add bowtie2-build-l bzip2recover glib-compile-resources mwm ldbdel grub-kbdcomp gifclrmp lstopo-no-graphics ipcmk gtbl pk-example-frobnicate coredumpctl htsfile filefrag nsec3hash host gdlib-config lzcat podchecker chpasswd qstatus vimdiff ip6tables-legacy sadf ip nslcd seal.sh peekfd audisp-remote cache_repair zstd xsubpp blkid wall ipcrm conf.py qrdel sge_shepherd splitnextera.sh filterqc.sh avahi-publish-service countgc.sh psicc idle3.7 python ifconfig rsync nl-route-get nsenter ldns-dpa tiffinfo smtp-sink sensors geqn uname cleanupdelta secret-tool ip6tables-save icupkg jfs_fsck db_tuner xsltproc nstat numactl ntptrace slappasswd sensors-conf-convert plotflowcell.sh dbus-test-tool mysql_config systemd-sysusers dhcpcd perlivp env_parallel.bash bbest.sh geoiplookup6 arptables-nft postlog less roff2html ntpq gobject-query rm tiffcp postmulti dedupe.sh prove nl-list-sockets migspeed ldns-signzone nfsiostat profiles dmsetup prtstat pdata_tools bc nl-link-set nfnl_osf gfind_missing_files dlist_test grpck dirmngr-client SHLIB w iusql make xfs_fsr pkcheck dnssec-checkds csplit samba-regedit Sweave smbcacls sha256sum unshare findssl.sh xzdiff egrep eject pod2html smbstatus Stangle as g++ gd2topng addftinfo setfont slirpvde cal expect bbmerge.sh sort soelim install-sh lzless pinky nl-route-list mount.nfs4 install mktemp pacat gdcmpgif era_dump hwloc-annotate troff samba_kcc swapon getopt thin_check tkpasswd cc gitk tiffsplit aws unlz4 gdtopng ndiff nf-exp-delete truncate mbni-dirfingerprint nameif tracepath iptables-legacy-restore jfs_mkfs pinentry-emacs avahi-dnsconfd bedpeToBam groupBy env fallocate lsinitcpio unicode2ascii.sh gtester estherfilter.sh icu-config tadpipe.sh msggrep taxserver.sh gpg-zip cut sgepasswd sge-disable-submits msguniq ldns-read-zone wirefilter pivot_root gradesam.sh setkeycodes bbduk.sh resizecons pldd clockdiff fusermount-glusterfs fsck.reiserfs lesskey hwloc-gather-topology xfs_freeze vde_tunctl grub-syslinux2cfg geoiplookup libassuan-config sdl2-config afmtodit shred hwloc-assembler msgattrib pam_tally2 gdparttopng showkey psfxtable pamon dir certtool switch_root base64 cmp invertkey.sh cd-iccdump bbcms.sh gluster-georep-sshkey ifcfg vde_switch fc-cat capsh xkibitz lnstat libtoolize dwebp removehuman2.sh partx cjpeg login rvim bash fmt jasper rnano Rscript tzselect colormgr ip6tables-restore-translate dmstats mcookie rst2pseudoxml.py umount.nfs encguess giftool javareconf lzdiff h2ph setmetamode gunzip find rsh webpinfo x86_64-pc-linux-gnu-gcc-ar smbclient rst2odt cutprimers.sh zipsplit pngtogd2 gifinto xfs_admin summarizemerge.sh shuf bamToFastq eventlogadm sserver telinit routef b2sum pstree.x11 xa2multi.pl showmount slapadd gpinyin avahi-discover-standalone mysql_plugin icuinfo randomgenome.sh asn1Coding df hpftodit multiIntersectBed mbni-sendthem calc_tickadj rcp rftp xzless gi2ancestors.sh python3 cutgff.sh hwloc-ps lsusb join zramctl sclient rst2html4 pzstd gluster usermod texi2pdf python3.7m hwloc-bind ndrdump git-cvsserver bashbug tclsh mdmon nscd sge_coshepherd modinfo rename.sh gnutls-cli gpg-error-config sync mev sleep nfsdcltrack env_parallel.zsh rst2xetex tapestat fdformat removemicrobes.sh talk head daisync-locate tsort build fc-list mknod filterbarcodes.sh ldapexop tdbtool ranlib lessecho chrt gdiffmk mkfs.reiserfs x86_64-pc-linux-gnu-gcc-8.3.0 nl-fib-lookup asn1Parser nl-neigh-list mpicalc thin_restore aureport pyvenv pvscan png2pnm INSTALL corelist kibitz netstat qrstat pyrsa-verify vte-2.91 gpgme-tool rev gpg-error mbni-ipcrm x86_64-pc-linux-gnu-g++ sar gxditview rst2html.py clusterBed dnssec-keymgr lslogins mkpasswd ldattach addr2line bedutil ldapsearch ldns-config fsadm unlink samtoroc.sh ksba-config flankBed k5srvutil ldbedit bed12ToBed6 ipmaddr fgconsole trietool-0.2 swaplabel nfs-cp stty unzipsfx dbus-monitor bssh unix2mac makepkg nettle-lfib-stream bzmore mkfs.ext4 ocsptool xtables-nft-multi zstdcat qping sed libhts.so qsub ldns-notify mysql_upgrade smbtorture prlimit chroot dmesg i686-pc-linux-gnu-pkg-config applygnupgdefaults gluster-setgfid2path tipc localectl dc ldns-compare-zones ctstat glib-compile-schemas pyrsa-sign uuidgen du rtmon mount.fuse cache_metadata_size grub-mknetdir xslt-config cifs.upcall luac nl-cls-add lsipc dos2unix mexext gpm ptardiff qemu-img size ebtables-nft-save postfix-collate.pl bsdcat ecpg qselect wrjpgcom era_invalidate shasum ssh-agent commonkmers.sh toe curl lvscan hb-shape nl-route-add postdrop xzgrep nl-neightbl-list automount filterlines.sh glurp rst2html setcap sortBed mysqlshow getfattr perldoc shiftBed sqlite3 modutil bootctl omxregister-bellagio instmodsh chown zstdless mysqlcheck gtk-query-settings zipdetails e2scrub lvm zsoelim db_verify xauth msgfilter linksBed dedupebymapping.sh mkfs.ext2 pcap-config annotateBed systemd-tmpfiles git-receive-pack pdbedit route lsmod aulast lvconvert vde_l3 lstopo dumpmscat sulogin pileup.sh krb5-config ausearch desktop-file-validate ownership bbwrap.sh psql update-desktop-database nl-link-enslave grn smbspool clusterdb umount sasldblistusers2 libhts.so.2 vgextend cxpm lsblk readelf getOverlap rst2latex delv keyctl drill rpc.nfsd lvmsar deallocvt BATCH db_load objdump db_checkpoint processhi-c.sh locale-gen ntpd hwloc-assembler-remote scmp_sys_resolver thin_metadata_size thin_ls gpgconf jp.py roff2pdf nl-addr-list pidl printf qshape fsck.ext4 namei testparm js52-config stat virtfs-proxy-helper hwloc-diff nfsconf LINK tr rarpd gpgtar sudo resizepart runcon lsns subsketch.sh pvremove cat getfacl gcron.py db_dump vwebp getent.ldap hwloc-ls sketchblacklist.sh screen-4.6.2 tabs tetramerfreq.sh avahi-publish groupdel veritysetup unzip wifi-menu addpart init chattr tset tiffcmp lzgrep agetty mac2unix Rdiff db_deadlock summarizequast.sh mergeBed slapd mysqltest comparevcf.sh env_parallel.csh reiserfsck halt sketch.sh nologin sge_shadowd sprof uil vde_plug2tap udevadm xmlsec1-config xtables-monitor pcre2grep cache_restore nl-neigh-add intersectBed x86_64-pc-linux-gnu-gcc-ranlib isosize pal2rgb demuxbyname2.sh lvremove chem tune2fs pkg-config column systemd-detect-virt nl-class-delete check nl-tctree-list systemctl sndfile-info cwebp blkdiscard djpeg grademerge.sh ldbmodify get_driver ebtables-nft rstpep2html fsck.ext2 bedToBam pkcs1-conv start-statd vimtutor genl-ctrl-list nucBed filterbycoverage.sh arp kproplog bbfakereads.sh taxtree.sh integritysetup write ldns-resolver systemd-ask-password filterbysequence.sh niceload unzstd gtk3-widget-factory callpeaks.sh ivshmem-client sum dnssec-cds mkfs asn1Decoding qresub bzgrep request-key grpconv unpigz isql kmerlimit2.sh kbd_mode pcre2-config whereis splitsam6way.sh roff2ps parec globusconnectpersonal kinit dbus-uuidgen showjournal signtool get_stat nmapfe smtp-source subtractBed mbni-run-vm unbuffer textfile.sh env_parallel.mksh ncat debugfs seq ldns-zcat cairo-sphinx xfs_spaceman cracklib-unpacker Rd2pdf audispd-zos-remote xfs_copy newgrp time gencmn awk khist.sh enc2xs oLschema2ldif vgdisplay env_parallel.pdksh ss objcopy iptables e2fsck cups-config grub-fstest env_parallel.dash mysqlimport xmlcatalog pinentry-curses bbmask.sh qemu-io rst2xml.py mkinitcpio su aws_zsh_completer.sh sshd mount.glusterfs summarizescafstats.sh wdctl aws_completer rst2xetex.py pkill bgzip R setleds runuser rpc.mountd mkfs.minix plotgc.sh jfs_logdump glusterfs cifs.idmap rarp bowtie2-build-l-debug nl-link-release pg_isready ip6tables-nft-save rqcfilter.sh machinectl get_device elf2dmp mklost+found dnssec-revoke repo-add repair.sh pcretest env_parallel.ksh xxd uname26 tiffgt visudo multiBamCov mysqlbinlog pactl bbrename.sh gtk-encode-symbolic-svg xfs_metadump iptables-legacy ivshmem-server hostid tiffmedian update-mime-database mkfs.bfs slapacl psl ldbadd blkmapd fsck.jfs pinentry python3.7m-config users xzmore iptables-restore hwloc-info dnssec-signzone pcre2test nf-ct-list texi2dvi qhost postmap padsp fsck.cramfs nm glusterfind rst2html4.py idle2 pmap watch postlock libnetcfg thin_dump vde_over_ns uuclient createuser infotocap gawk znew sndfile-metadata-set parecord ip6tables-nft-restore vdecmd mergesorted.sh rdjpgcom sndfile-cmp test config.status removehuman.sh rmmod chmem pvcreate rtcwake comparesketch.sh git-upload-archive sortbyname.sh hltest jemalloc.sh tiffset xfs_io cpp ip6tables-legacy-restore readprofile pvs iptables-xml pango-list xml2-config activate_matlab.sh fsck pgrep kprop lzmadec cd-create-profile grap2graph version.sh pod2usage iostat grub-sparc64-setup pydoc2 sysctl posttls-finger process-scheduler-log nslookup grub-glue-efi rpc.idmapd nroff comm db_hotbackup nfsstat filterbyname.sh qualfa2fq.pl synthmda.sh c89 idle-nodes xfs_growfs derb roff2text idtree.sh realpath unixcmd mkreiserfs sntp uuidparse pairToPair dbiprof xnmap strace-log-merge 2to3-2.7 losetup dbiproxy f95 jeprof slattach ldns-version xfs_mkfile passwd aserver chacl pvchange linkicc rlogind qstat strings mountpoint loadkeys parallel sql pwmconfig ldns-chaos gtester-report vdekvm testformat2.sh ldapdelete dnsdomainname pydoc3.7 gifecho qtcsh segment resolvconf nice qmod ddns-confgen dmeventd ldapcompare jfs_tune fsck.minix tknewsbiff iptables-nft-restore msgfmt lvs sendsketch.sh mandb removebadbarcodes.sh vmstat manpath samba_upgradedns linux32 hwloc-distances mergebarcodes.sh taskset iptables-save autoexpect lvcreate ionice desktop-file-install pwd imginfo edit rlogin mkfs.jfs bshell podselect datamash tiff2bw removecatdogmousehuman.sh sudoedit memhog rtacct lvchange mbni-reserve-nodes mkfs.cramfs xgettext lzmore lvmdump guile systemd-mount callvariants2.sh diff3 regpatch bunzip2 gdbm_dump processspeed.sh wipefs xmllint get-versions smbcontrol desktop-file-edit cdb lkbib pic rst2s5.py splitsam4way.sh libpng16-config grub-render-label chcpu avahi-autoipd fc-match ktutil bowtie2-build-s-debug jiv tiffdither ex cifscreds split ftpd lslocks gensprep parset rst2odt_prepstyles.py psfaddtable bowtie2-align-l texi2any bvnc isadump tickadj key.dns_resolver talkd xzcmp tiff2rgba badblocks bwa systemd-nspawn slopBed avahi-resolve-address fstrim slabtop isaset iptables-legacy-save qsched vgcfgrestore roff2dvi db_recover perl nl-cls-list autopasswd samba-tool ldapadd cache_check lua matlab-glselector.sh sndfile-concat nl-pktloc-lookup ldns-rrsig vipw addgnupghome transicc vi mkhomedir_helper gr2fonttest libtool cvtsudoers auditctl zstdmt lvmetad cryptsetup crossblock.sh vde_plug rst2latex.py explodetree.sh parcat dbus-launch webpng vgsplit certutil bbnorm.sh msgcat REMOVE gpgme-config sge-enable-submits grolj4 systemd-inhibit renice msgcomm e2freefrag pkexec ldns-gen-zone iptables-translate xfs_repair printenv rename mysqldump gdk-pixbuf-csource vdeq newgidmap makepolymers.sh hmac256 giffix vgrename bbmap.sh resolvectl postalias avahi-daemon perlbug getconf mesg uniq gitable.sh mtrace bzdiff groupmod piconv xfs_bmap vim bedToIgv nfbpf_compile chgrp e2undo gcov-tool qmqp-source passmass dnssec-importkey funzip Rcmd ldns-key2ds rst2man.py lzma vgreduce ldns-mx fungalrelease.sh zmore grub-mkstandalone setterm vgck ntptime look img2webp testpkg python2 sharesec basenc groups decryptdir setcifsacl zdump smbcquotas gifsponge bbrealign.sh who trietool shuffle2.sh systemd-stdio-bridge ldconfig ip6tables-legacy-save grpunconv gcc-ar lzcmp pacman-key lvmsadc kbdrate msgexec mapscrn tftpd shuffleBed gpg-agent sndfile-interleave named-rrchecker globusconnect jfs_debugfs thin_repair logsave ldappasswd mdb_dump nmbd bowtie2-align-s Rdconv qemu-system-x86_64 jobstats dbus-update-activation-environment userdel gpm-root ldns-update sem systemd-resolve xpstat systemd-machine-id-setup false bzip2 tunefs.reiserfs ps lvmconfig nl krb5-send-pr pdftexi2dvi wbinfo idle install-info pyrsa-priv2pub unix_chkpwd tadpole.sh pngtogd gtk-update-icon-cache e2image nl-qdisc-delete gi2taxid.sh zgrep grub-reboot file grub-ofpathname whoami gzexe nf-monitor flac mdb_stat newusers zic blockdev mbni-receiveit dnssec-keyfromlabel reformat.sh taxonomy.sh pcre-config eqn xfs_mdrestore paste pwconv escapesrc busctl daisync cache_writeback splitsam.sh useradd ssh-add rpcinfo ip6tables-restore ldapwhoami auditd vgimportclone sge_share_mon named-journalprint partition.sh ls tsig-keygen pinentry-qt preconv groff sxpm tcpdump rpcclient weather vlock systemd-cat symkeyutil mailx demuxbyname.sh kdb5_util configure tic lookbib readlink col showdb fgrep regdiff touch dirname xtrace reboot cfdisk vedit bbsplit.sh pscap tiff2ps postqueue gnutls-cli-debug consect.sh ld.bfd bdf2gdfont.pl x86_64-pc-linux-gnu-gfortran which chfn memusagestat smbget nl-qdisc-add cracklib-format ldnsd bzcat guile-config vgremove i386 mmroff gifbuild reiserfstune summarizesketch.sh pam_tally apropos indxbib findsmb plurp findfs tree pkgfile portablectl sndfile-salvage ld.gold npth-config createdb gencnval gettext guile-snarf smbpasswd depmod expr splain wget pacman-db-upgrade numastat systemd-path pwunconv vigr pkgdata fsck.xfs grub-probe mutate.sh gsettings xfs_info tdbrestore delpart sndfile-play mapPacBio.sh tcsh systemd-id128 matrixtocolumns.sh fold openvt display-buttons bloomfilter.sh free timeout ldns-keyfetcher genbrk tadwrapper.sh hwloc-calc ncursesw6-config gc.py qemu-nbd canberra-boot avahi-discover rdisc xfs_scrub_all scriptreplay hwloc-distrib lvextend yes slapschema signver uconv fsck.ext3 infocmp mountstats ipcs ldbrename xzegrep sim_client gperl giftogd2 pfbtops grub-mkfont systemd-escape poweroff showstat4 cmsutil tdbdump odbc_config gropdf kmerfilterset.sh nl-rule-list vgconvert dead-nodes nping db_upgrade lastlog wc update-leap cryptsetup-reencrypt sm-notify strip nfsidmap auvirt libpng-config smartd zipgrep glustereventsd accessdb tabix gcc summarizecrossblock.sh mount.nfs avahi-resolve-host-name gendict expiry update-pciids gpgscm python3-config lpunlock ldns-keygen hwloc-compress-dir utmpdump findmnt systemd-hwdb odbcinst blkdeactivate kmercoverage.sh pinentry-gtk-2 gtk-query-immodules-3.0 STARlong renameimg.sh metaflac tput zip readqc.sh nano db_log_verify python2.7-config slapdn vercmp raw ulockmgr_server unix2dos vgscan unexpand nl-util-addr gzip resize2fs qconf avahi-bookmarks tdbbackup png-fix-itxt unicode_stop smartctl rst2html5.py stdbuf grub-script-check last translate6frames.sh kill slapindex ssh-keyscan qevent getkeycodes logrotate hwloc-dump-hwdata grub-mklayout bsdcpio vgimport nmblookup tclsh8.6 mkfs.xfs exportfs arpd bdftogd date removesmartbell.sh lz4 pstree pidof mpstat makedb mysqlslap netcap fc-query pinentry-gnome3 2to3 gennorm2 mariadb_config sha1sum regtree ar ctrlaltdel dumpkeys ntp-wait vcf2gff.sh dislocate 2to3-3.7 ldns-dane cracklib-check usbhid-dump grub-set-default gtk3-icon-browser masktest lzegrep post-grohtml giftext gpgv2 uuserver mkfs.ext3 srptool gresource pydoc3 expand idmatrix.sh kpropd roff2x gio-querymodules webcheck.sh grolbp env_parallel.sh rsvg-convert event_rpcgen.py pkg-config-32 luac5.3 slapcat xfs_ncheck systemd-cgtop systemd-delta rtstat chgpasswd dirmngr tload bowtie2 multixterm mapBed bbcountunique.sh proxy unix_update sensord x86_64-pc-linux-gnu-gcc calcmem.sh rclone tfmtodit vdir nl-qdisc-list qalter summarizeseal.sh rtpr env_parallel.ash tee annotate bowtie2-build bedtools pnm2png celtdec051 loadreads.sh recode-sr-latin glib-gettextize samba_spnupdate dbus-send tificc qrsub expandCols callvariants.sh gio csh xz base32 sh captest gprof vde_cryptcab gpgv catchsegv yat2m vacuumdb uuidd qhold nohup compile_et psktool hb-view rpc.gssd README.md pam_timestamp_check zforce dbus-cleanup-sockets tbl fc-scan setvtrgb [ iptables-restore-translate gusbcmd mount kdb5_ldap_util mex glusterfsd makecontaminatedgenomes.sh pango-view lua5.3 rtags db_stat decontaminate.sh chsh.ldap logger ntp-keygen localedef dedupe2.sh zdiff tc cifsdd reset lcdata.xsd glusterd mergeOTUs.sh rdma pidstat openssl hb-subset timedatectl grub-mkconfig nproc rstpep2html.py cracklib-packer oathtool xzdec glfsheal cache_dump ldns-revoke ip6tables-nft dircolors blkzone wayland-scanner xtables-legacy-multi slaptest strace grub-editenv strace-graph windowBed lsattr qselect-node-list systemd-tty-ask-password-agent makeinfo dbus-daemon lvmconf kbdinfo gentest pvmove msa.sh gpg-connect-agent qmake gtk3-demo gpg-wks-server nfs-cat usb-devices raw2tiff ethtool grops dropdb shutdown echo e2label e2mmpstatus era_check vdeqemu msginit maskFastaFromBed git-upload-pack genccode python2.7 lastb ptx representative.sh nodes-in-job setfacl license.txt regshell gif2rgb sln krb5kdc pcregrep bbmerge-auto.sh gpg cd-it8 bowtie2-inspect-s msgunfmt mkfifo lexgrog gss-server systemd-run kpasswd croco-0.6-config memusage gif2webp rshd groffer git fincore ldapurl isc-config.sh kcompress.sh chage nf-exp-list filterbytile.sh env_parallel.fish slurp avahi-publish-address xmlsec1 getcap msgconv cd-fix-profile attr fax2tiff matlab setfattr x86_64-pc-linux-gnu-pkg-config jpegtran rst2odt.py ebtables-nft-restore gpgme-json bowtie2-inspect-s-debug dnssec-dsfromkey lsof vgmknodes repo-elephant statswrapper.sh grub-mkrelpath virgl_test_server c99 gpgparsemail kmod ld printtime.sh tjbench newaliases crosscontaminate.sh python-config gapplication ppm2tiff xfs_scrub trust samba-gpupdate kvno nl-list-caches ldns-test-edns ftp-rfc randomBed mergesam.sh xmlwf windowMaker zless named-nzd2nzf xfs_db captoinfo gdbmtool dpipe stats.sh postfix postfilter.sh p11-kit qlogin db_replicate pydoc top nettle-pbkdf2 gpgsm kswitch canberra-gtk-play samtools randomreads.sh systool systemd-firstboot chmod debugfs.reiserfs argon2 dnssec-coverage ssh-copy-id timed-read avahi-browse unionBedGraphs biosdecode unxz pod2texi dnssec-settime clumpify.sh mdb_copy kmercountexact.sh mk_cmds qacct jfs_fscklog sdiff groupadd symcryptrun unlzma nl-link-list grub-mkpasswd-pbkdf2 routel systemd-analyze dbwrap_tool idiag-socket-details gtk3-demo-application fsfreeze Rprof kapastats.sh numfmt fuse.sh getent);
$| = 1;
require "$qDir/lib/main.pl";
qMain();
