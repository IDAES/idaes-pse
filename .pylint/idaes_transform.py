"""
Pylint transform plugin to make Pylint aware of custom `ProcessBlock` classes
dynamically created through the `declare_process_block_class()` decorator.

See #1159 for more information.
"""
import functools
import logging

import astroid
from astroid.builder import extract_node, parse


_logger = logging.getLogger('pylint.ideas_plugin')


# TODO figure out a better way to integrate this with pylint logging and/or verbosity settings
_display = _notify = print


def has_declare_block_class_decorator(cls_node, decorator_name="declare_process_block_class"):
    if 'idaes' not in cls_node.root().name:
        return False
    decorators = cls_node.decorators
    if not decorators:
        return False
    for dec_subnode in decorators.nodes:
        if hasattr(dec_subnode, 'func'):
            # this is true for decorators with arguments
            return dec_subnode.func.as_string() == decorator_name
    return False


# this will be called N times if running with N processes
# returning the cached result on subsequent calls
@functools.lru_cache(maxsize=1)
def get_base_class_node():
    _notify('Getting base class node')
    import_node = extract_node('from idaes.core.process_block import ProcessBlock; ProcessBlock')
    cls_node = next(import_node.infer())
    return cls_node


def add_attribute_nodes(node: astroid.ClassDef, attr_names):
    for attr_name in attr_names:
        rhs_node = astroid.Unknown(
            lineno=node.lineno,
            parent=node,
        )
        node.locals[attr_name] = [rhs_node]
        node.instance_attrs[attr_name] = [rhs_node]


def create_declared_class_node(decorated_cls_node: astroid.ClassDef):
    decorators = decorated_cls_node.decorators
    call = decorators.nodes[0]
    name_arg_node = call.args[0]
    decl_class_name = name_arg_node.value

    base_class_node = get_base_class_node()
    decl_class_node = astroid.ClassDef(
        decl_class_name,
        # TODO the real doc should be available as the "doc" kwarg of the decorator
        # but it's not clear if we're going to need it anyway
        doc=f"Declared from {decorated_cls_node.name}",
        parent=decorated_cls_node.parent,
        lineno=decorated_cls_node.lineno,
    )
    decl_class_node.bases.extend(
        [
            base_class_node,
            decorated_cls_node,
        ]
    )
    # doesn't seem to be needed at the moment
    # add_attribute_nodes(decl_class_node, ['_orig_name', '_orig_module'])
    return decl_class_node


def is_idaes_module(mod_node: astroid.Module):
    mod_name = mod_node.name
    if mod_name:
        _display(f'analyzing module: {mod_name}')
    return 'idaes' in mod_name


def register_process_block_class(decorated_cls_node: astroid.ClassDef):
    module_node = decorated_cls_node.parent
    _display(module_node.name)
    decl_class_node = create_declared_class_node(decorated_cls_node)


def iter_process_block_data_classes(mod_node: astroid.Module):
    for node in mod_node.body:
        if isinstance(node, astroid.ClassDef):
            if has_declare_block_class_decorator(node):
                yield node


def register_process_block_classes(mod_node: astroid.Module):
    _display(mod_node.name)
    for decorated_cls_node in iter_process_block_data_classes(mod_node):
        decl_cls_node = create_declared_class_node(decorated_cls_node)
        _display(f'{decorated_cls_node.name} -> {decl_cls_node.name}')


def is_base_pyomo_var_class(node):
    try:
        _display(f'node.qname()={node.qname()}')
        return 'pyomo.core.base.var.Var' in node.qname()
    except AttributeError:
        pass
    return False


def get_concrete_pyomo_var_class(node: astroid.ClassDef, context=None) -> astroid.ClassDef:
    # node_from_import = extract_node('from pyomo.core.base.var import SimpleVar; SimpleVar')
    # simple_var_cls_node = astroid.helpers.safe_infer(node_from_import)
    clsdef_code = """
    from pyomo.base.var import SimpleVar, IndexedVar
    class ConcreteVar(SimpleVar, IndexedVar):
        pass
"""
    simple_var_cls_node = extract_node(clsdef_code)
    print(simple_var_cls_node)
    return simple_var_cls_node


def infer_concrete_var_instance(node: astroid.ClassDef, context=None):
    _display(f'abstract var class: {node}')
    concrete_cls_node = get_concrete_pyomo_var_class(node, context=context)
    _display(f'concrete var class: {concrete_cls_node}')
    return iter([concrete_cls_node.instantiate_class()])


def register(linter):
    "This function needs to be defined for the plugin to be picked up by Pylint"


astroid.MANAGER.register_transform(
    # NOTE both these options were tried to see if there was a difference in performance,
    # but it doesn't seem to be the case at this point
    # astroid.ClassDef, register_process_block_class, has_declare_block_class_decorator
    astroid.Module, register_process_block_classes, is_idaes_module,
)

astroid.MANAGER.register_transform(
    astroid.ClassDef, astroid.inference_tip(infer_concrete_var_instance), is_base_pyomo_var_class
)
