Running Mkbootstrap for IndexedBase_168b ()
chmod 644 "IndexedBase_168b.bs"
"/usr/bin/perl" -MExtUtils::Command::MM -e 'cp_nonempty' -- IndexedBase_168b.bs blib/arch/auto/Bio/DB/IndexedBase_168b/IndexedBase_168b.bs 644
"/usr/bin/perl" "/usr/share/perl5/ExtUtils/xsubpp"  -typemap "/usr/share/perl5/ExtUtils/typemap"   IndexedBase_168b.xs > IndexedBase_168b.xsc
mv IndexedBase_168b.xsc IndexedBase_168b.c
gcc -c  -I"/opt/scripts" -D_REENTRANT -D_GNU_SOURCE -fno-strict-aliasing -pipe -fstack-protector -I/usr/local/include -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -O2 -g -pipe -Wall -Wp,-D_FORTIFY_SOURCE=2 -fexceptions -fstack-protector --param=ssp-buffer-size=4 -m64 -mtune=generic   -DVERSION=\"0.00\" -DXS_VERSION=\"0.00\" -fPIC "-I/usr/lib64/perl5/CORE"   IndexedBase_168b.c
In file included from /usr/include/sys/types.h:30:0,
                 from /usr/lib64/perl5/CORE/perl.h:605,
                 from IndexedBase_168b.xs:2:
/usr/include/bits/types.h:179:12: error: unknown type name ‘__FSWORD_T_TYPE’
 __STD_TYPE __FSWORD_T_TYPE __fsword_t;
            ^
/usr/include/bits/types.h:184:12: error: unknown type name ‘__SYSCALL_SLONG_TYPE’
 __STD_TYPE __SYSCALL_SLONG_TYPE __syscall_slong_t;
            ^
/usr/include/bits/types.h:186:12: error: unknown type name ‘__SYSCALL_ULONG_TYPE’
 __STD_TYPE __SYSCALL_ULONG_TYPE __syscall_ulong_t;
            ^
IndexedBase_168b.xs: In function ‘_strip_crnl’:
IndexedBase_168b.xs:11:27: warning: value computed is not used [-Wunused-value]
         for (s = str; *s; *s++) {
                           ^
make: *** [IndexedBase_168b.o] Error 1
