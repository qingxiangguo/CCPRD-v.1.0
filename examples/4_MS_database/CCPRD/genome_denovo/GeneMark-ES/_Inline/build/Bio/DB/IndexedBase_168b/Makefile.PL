use ExtUtils::MakeMaker;
my %options = %{
{
  'TYPEMAPS' => [
    '/usr/share/perl5/ExtUtils/typemap'
  ],
  'NAME' => 'Bio::DB::IndexedBase_168b',
  'INC' => '-I"/opt/scripts"',
  'VERSION' => '0.00'
}
};
WriteMakefile(%options);

# Remove the Makefile dependency. Causes problems on a few systems.
sub MY::makefile { '' }
