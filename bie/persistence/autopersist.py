#!/usr/bin/env python

import sys, os, re, string
from string import Template

# legal C++ identifier (regexp)
c_id = '[\s]*\w*[\s]*'

# single word identifier. i.e. "int", "long", "hello_world", etc.
# groups the identifier for later retrieval
id = '[\s]*(\w*)[\s]*'

# multi word identifier. i.e. "unsigned int".
# groups the whole identifier as one for later retrieval
multi_id = '[\s]*([\w ]*)[\s]*'

all_classes = []
classes = []
header = ''

suffix = '.ap'

serial_base_template = Template(
    '$indent try {                                                         \n' +
    '$indent  ar $op BOOST_SERIALIZATION_BASE_OBJECT_NVP($var);            \n' +
    '$indent  BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     \n' +
    '$indent }                                                             \n')

serial_template = Template(
    '$indent try {                                                         \n' +
    '$indent  ar $op BOOST_SERIALIZATION_NVP($var);                        \n' +
    '$indent  BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     \n' +
    '$indent }                                                             \n')


class persistent_class:
    """
    Using reflection and regexp, annotations that start with @ are
    matched automatically by the corresponding at_<name>
    function. e.g. upon matching @persistent(Node), at_persistent()
    function will be invoked and the string "Node" will be passed in
    as the second (after self) parameter.

    For exceptional cases, like namespace matching, the mapping of the
    match regexp and the name of the function needs to be specified
    manually in the "custom_match"" attribute.

    It is important to note that the strings found in parentheses
    within the regexp are picked up and passed in as parameters to the
    corresponding function. It is _your_ responsibility to match the
    number of such subexpressions found in regexp to the number of
    formal parameters of the corresponding function. Failing to do so
    won't raise exception. It just won't match and it will be in the
    output, unaltered. Ordering is preserved.
    
    For example, you'd probably want to enclose the entire regexp in a
    pair of parenthesis if it is to be matched with "verbatim"
    function, as verbatim function simply emits what gets passed in,
    unchanged.
    """
    includes = [r'"Serializable.h"']
    custom_match = \
        {
        # fast-trap namespace keyword
        '(using namespace)': 'verbatim',
        'namespace ' + id: 'process_namespace',
        '@IMPLICIT_DEFAULT_CONSTRUCTOR()':
            'process_default_constructor',
        # match pure virtual functions. First matches "virtual", then
        # everything except ; or {. Once it finds { before ;, it bails
        # out since it's clearly not a pure virtual function. Once it
        # finds ;, perform a lookbehind match for 0. If 0 is the last
        # char right before ;, then this is a pure virtual function.
        '([ \t]*virtual[^\;\{]*(?<=0)\s*\;)' : 'process_pure_virtual_function',
        # handling two types of comments
        '(\/\/.*\n)': 'verbatim',
        '(\/\*(?:.|[\n])*?\*\/)': 'verbatim',
        # match the end of the file --> BOOST_CLASS_TYPE_INFO
        '(\#endif\s*$$)': 'generate_type_info',
        # wrap BOOST_SERIALIZATION lines
        '(.*BOOST_SERIALIZATION.*\[^\\\]\n)': 'wrap_boost_serialization'}
    marker = \
        '// AUTO GENERATED BY ' + sys.argv[0] + '\n'
    persistence_end_intro = \
        [ 'private:\n',
          'friend class boost::serialization::access;\n']
    serialize_sig = \
        [ 'template<class Archive>\n',
          'void serialize(Archive & ar, const unsigned int version) {\n' ]
    save_sig = \
        [ 'template<class Archive>\n',
          'void save(Archive & ar, const unsigned int version) const {\n' ]
    load_sig = \
        [ 'template<class Archive>\n',
          'void load(Archive & ar, const unsigned int version) {\n' ]

    def __str__(self):
        return ' | '.join([self.namespace, ' '.join(self.super), self.classname])

    def init_classvar(self):
        self.classname = ''
        self.super, self.regular_field, self.array_field = [], [], []
        self.def_ctor, self.is_abstract, self.is_base = False, False, True
        self.indent = ''

    def compile_match(self):
        self.match_map.update(self.custom_match)
    #  process match_map to produce a dictionary of..
    # { <compiled regexp pattern obj> : <bound function object>, ... }
        self.match_map_compiled = \
            dict([(re.compile(key), getattr(self, self.match_map[key])) 
                  for key in self.match_map.keys()])

    def __init__(self, fromFile):
        self.fromFile = open(fromFile, 'r')
        self.toFile = open(fromFile.replace(suffix, ''), 'w')
        self.namespace = ''
        self.init_classvar()
        # generate appropriate regexp for methods that start with at_
        def gen_regex(func_name):
            ret = func_name.replace('at_', '(\@')
            argcount = getattr(self, func_name).func_code.co_argcount
            if (argcount-2):
                ret += '\('
                for i in range(argcount-3):
                    ret += multi_id + ','
                ret += multi_id + '\)'
            ret += ')'
            return ret
        # reflectively search for at_ functions and create appropriate
        # regex to match them
        self.match_map = \
            dict([(gen_regex(x),x) for x in 
                  filter(lambda x: x.find('at_') == 0, dir(self))])
        # for debugging, uncomment the following to see the generated regex
        # print self.match_map.keys()
        self.compile_match()

# for convenience
    
    def out_plain(self, str):
        self.toFile.write(str)

    def out(self, str):
        if (re.compile('^[ \t]+').search(str)
            or len(str) == 0):
            self.out_plain(str)
        else:
            self.out_plain(self.indent + str)
            
    def outlist(self, list):
        for line in list:
            self.out(line)

    def out_indent(self):
        self.indent += '    '

    def out_unindent(self):
        self.indent = self.indent[4:]

    def print_match_map(self):
        mmap = self.match_map_compiled
        print '\n'.join(map(lambda x:x.pattern+' -> '+mmap[x].im_func.func_name,
                            mmap.keys()))
        
    def set_classname(self, name):
        self.classname = name
        # add a rule for detecting default constructor
        # 1. match destructor first to get it out of the way
        self.custom_match.update(
            {'(~'+name+'[ ]*\([ ]*\)[ ]*)':'verbatim'})
        # 2. match default constructor
        self.custom_match.update(
            {'([ ]'+name+'[ ]*\([ ]*\)[ ]*)':'process_default_constructor'})
        self.compile_match()
        
# business logic goes in here - what to do when there's a match
    def at_include_persistence(self, orig):
        for include in self.includes:
            self.out("#include " + include + "\n")

    def at_persistent_root(self, orig, name):
        self.set_classname(name)
        self.super.append('Serializable')
        self.out_plain(name + ': public Serializable')

    def at_persistent(self, orig, name):
        self.set_classname(name)
        self.out_plain(name)

    def at_super(self, orig, name):
        self.is_base = False
        self.super.append(name)
        self.out_plain(name)

    def at_autopersist(self, orig, name):
        self.regular_field.append(name)
        self.out_plain(name)

    def at_autopersist_array(self, orig, type, name, size):
        self.array_field.append((type, name, size))

    def at_no_autopersist(self, orig, name):
        self.out_plain(name)
        
    # abbreviated versions
    def at_ap(self, orig, name):
        self.at_autopersist(orig, name)

    def at_ap_array(self, orig, type, name, size):
        self.at_autopersist_array(orig, type, name, size)

    def at_no_ap(self, orig, name):
        self.out_plain(name)

    def persistent_end_prologue(self, orig):
        self.out(self.marker)
        if not self.def_ctor:
            self.out('protected:\n')
            self.out(self.classname+'() {}\n')
        self.outlist(self.persistence_end_intro)

    def persistent_end_epilogue(self, orig):
        # Export all classes
        if self.namespace:
            all_classes.append(self.namespace + '::' + self.classname)
        else:
            all_classes.append(self.classname)
        # Register this class only if it's not purely virtual or a base class
        if not self.is_abstract:
            if self.namespace:
                classes.append(self.namespace + '::' + self.classname)
            else:
                classes.append(self.classname)

        # initialize class variables since we're done with this class
        self.init_classvar()
        
    def at_persistent_end_split(self, orig):
        self.out("#ifndef SWIG\n")
        self.out_indent()
        self.persistent_end_prologue(orig)
        self.generate_split_member(orig)
        self.out_unindent()
        self.persistent_end_epilogue(orig)
        self.out_indent()
        self.out("#endif\n")
        self.out_unindent()
        
    def at_persistent_end(self, orig):
        self.out("#ifndef SWIG\n")
        self.out_indent()
        self.persistent_end_prologue(orig)
        # Determine if we need to split serialize()
        if len(self.array_field):
            # need to split if we have array_fields
            self.generate_split_member(orig)
        else:
            self.generate_serialize(orig)
        self.out_unindent()
        self.persistent_end_epilogue(orig)
        self.out_indent()
        self.out("#endif\n")
        self.out_unindent()
        
    def generate_serialize(self, orig):
        self.outlist(self.serialize_sig)
        self.out_indent()
        self.out('this->pre_serialize(ar, version);\n')
        # print BOOST_SERIALIZATION_BASE_OBJECT_NVP
        map(lambda x: self.out(serial_base_template.
                               substitute(indent=self.indent,
                                          op='&', var=x)), self.super)
        # print BOOST_SERIALIZATION_NVP
        map(lambda x: self.out(serial_template.substitute(indent=self.indent,
                                                          op='&',
                                                          var=x)),
            self.regular_field)
        self.out('this->post_serialize(ar, version);\n')
        self.out_unindent()
        self.out('}\n')

    def generate_split_member(self, orig):
        self.out('BOOST_SERIALIZATION_SPLIT_MEMBER();\n\n')
        self.generate_save(orig)
        self.generate_load(orig)

    def generate_catch_exception_save_load(self):
        self.out_unindent()
        self.out('} catch (::boost::archive::archive_exception & e) {\n')
        self.out_indent()
        self.out('throw new BoostSerializationException(string(e.what()),')
        self.out_plain('__buf, __FILE__,__LINE__);\n')
        self.out_unindent()
        self.out('} catch (BoostSerializationException * e) {\n')
        self.out_indent()
        self.out('throw new BoostSerializationException(e, __buf,__FILE__,')
        self.out_plain('__LINE__);\n')
        
    def generate_save(self, orig):
        self.outlist(self.save_sig)
        self.out_indent()
        self.out('this->pre_save(ar, version);\n')
        # print BOOST_SERIALIZATION_BASE_OBJECT_NVP
        map(lambda x: self.out(serial_base_template.
                               substitute(indent=self.indent,
                                          op='<<', var=x)), self.super)
        # print BOOST_SERIALIZATION_NVP
        map(lambda x: self.out(serial_template.substitute(indent=self.indent,
                                                          op='<<',
                                                          var=x)),
            self.regular_field)
        if len(self.array_field):
            self.out('char __buf[128];\n')
        for var in self.array_field:
            self.out('for(uint32 __i=0; __i<'+var[2]+'; __i++) {\n')
            self.out_indent()
            self.out('sprintf(__buf, \"'+var[1]+'_%d\", __i);\n')
            self.out('try {\n')
            self.out_indent()
            self.out('ar <<  boost::serialization::make_nvp(__buf, '+var[1]+'[__i]);\n')
            self.generate_catch_exception_save_load()
            self.out_unindent()
            self.out('}\n')
            self.out_unindent()
            self.out('}\n')
        self.out('this->post_save(ar, version);\n')
        self.out_unindent()
        self.out('}\n\n')
    
    def generate_load(self, orig):
        self.outlist(self.load_sig)
        self.out_indent()
        self.out('this->pre_load(ar, version);\n')
        # print BOOST_SERIALIZATION_BASE_OBJECT_NVP
        map(lambda x: self.out(serial_base_template.
                               substitute(indent=self.indent,
                                          op='>>', var=x)), self.super)
        # print BOOST_SERIALIZATION_NVP
        map(lambda x: self.out(serial_template.substitute(indent=self.indent,
                                                          op='>>',
                                                          var=x)),
            self.regular_field)
        if len(self.array_field):
            self.out('char __buf[128];\n')
        for var in self.array_field:
            self.out(var[1]+' = new '+var[0]+'['+var[2]+'];\n')
            self.out('for(uint32 __i=0; __i<'+var[2]+'; __i++) {\n')
            self.out_indent()
            self.out('sprintf(__buf, \"'+var[1]+'_%d\", __i);\n')
            self.out('try {\n')
            self.out_indent()
            self.out('ar >> boost::serialization::make_nvp(__buf, '+var[1]+'[__i]);\n')
            self.generate_catch_exception_save_load()
            self.out_unindent()
            self.out('}\n')
            self.out_unindent()
            self.out('}\n')
        self.out('this->post_load(ar, version);\n')
        self.out_unindent()
        self.out('}\n')
        
    def process_namespace(self, name):
        self.namespace = name
        self.out('namespace '+name+' ')

    def process_default_constructor(self, line):
        self.def_ctor = True
        self.out(line)

    def process_pure_virtual_function(self, line):
        # print "PURE VIRTUAL FUNCTION DETECTED"
        # print line
        self.is_abstract = True
        self.out(line)

    def verbatim(self, orig):
        self.out(orig)

    def generate_type_info(self, orig):
        self.out("#ifndef SWIG\n")
        map(lambda x: self.out('BIE_CLASS_TYPE_INFO('+x+')\n'), classes)
        # For BOOST >= 1.42
        map(lambda x: self.out('BIE_CLASS_EXPORT_KEY('+x+')\n'), all_classes)
        self.out("#endif\n")
        self.out(orig)
        
    def wrap_boost_serialization(self, orig):
        self.out("#ifndef SWIG\n")
        self.out(orig)
        self.out("#endif\n")

    # The heart of matching logic. Lines in the input file are
    # concatenated to form a giant string. Then regexp/s are matched
    # from the beginning.
    def process(self):
        self.inside_comment = False
        scan = ''.join(self.fromFile)
        while scan:
            matches = []
            for regex in self.match_map_compiled.keys():
                match = regex.search(scan)
                if match:
                    matches.append(match)
            # one that matches first. Ties broken by the length (longest)
            def comparator(x, y):
                if x.start() < y.start():
                    return -1
                elif x.start() > y.start():
                    return 1
                else:
                    if x.end() > y.end():
                        return -1
                    else:
                        return 1
            matches.sort(comparator)
            
            if not matches:
                # no more matches
                self.out_plain(scan)
                return
            else:
                m = matches[0]
                self.out_plain(scan[:m.start()])
                # call the corresponding function
                self.match_map_compiled[m.re](*list(m.groups()))
                scan = scan[m.end():]
        pc.toFile.flush()
        
# python style main function. This is what gets called if this script
# is called as a standalone, i.e. not as a part of something else.
# check persistence/Makefile.autopersist.am
if __name__ == '__main__':
    toProcessHeader = sys.argv[1]
    header = sys.argv[1].replace(suffix,'')
    pc = persistent_class(toProcessHeader)
    pc.process()

