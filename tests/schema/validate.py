"""validate.py: an in-process JSON Schema check for three Python test files.

The check CI and the e2e cases run is `og dev schema SCHEMA DOC`
(rust/crates/og-dev/src/schema.rs). This module stays only because
tests/test_collect_summary.py, tests/test_hla_report.py and
tests/test_render_report.py call validate(schema, doc) in-process, in the
unit-tests job, which has no og. It goes when those tests move to Rust.

It follows og dev schema's rules: the keywords tests/schema/summary.schema.json
uses (type, enum, required, properties, additionalProperties as a schema, true
or false, items, minItems, maxItems, minimum, allOf and local $ref); a bool is
never an integer or a number, also in enum; NaN fails minimum. An unknown
keyword, type name or $ref target raises SchemaError where the document
reaches it (og dev schema checks the whole schema before it reads the
document), so a check never passes because a rule was silently ignored.
"""
import math

KNOWN = {"$schema", "title", "description", "definitions", "type", "enum", "required", "properties",
         "additionalProperties", "items", "minItems", "maxItems", "minimum", "allOf", "$ref"}
TYPES = {
    "object": lambda v: isinstance(v, dict),
    "array": lambda v: isinstance(v, list),
    "string": lambda v: isinstance(v, str),
    "integer": lambda v: isinstance(v, int) and not isinstance(v, bool),
    "number": lambda v: isinstance(v, (int, float)) and not isinstance(v, bool),
    "boolean": lambda v: isinstance(v, bool),
    "null": lambda v: v is None,
}


class SchemaError(Exception):
    pass


def resolve(root, ref):
    if not ref.startswith("#/"):
        raise SchemaError(f"only local $ref is supported: {ref}")
    node = root
    for part in ref[2:].split("/"):
        if not isinstance(node, dict) or part not in node:
            raise SchemaError(f"$ref {ref} points at nothing")
        node = node[part]
    return node


def same(a, b):
    """Equal by JSON type: a bool never equals a number."""
    if isinstance(a, bool) or isinstance(b, bool):
        return isinstance(a, bool) and isinstance(b, bool) and a == b
    if isinstance(a, list) and isinstance(b, list):
        return len(a) == len(b) and all(same(x, y) for x, y in zip(a, b))
    if isinstance(a, dict) and isinstance(b, dict):
        return a.keys() == b.keys() and all(same(a[k], b[k]) for k in a)
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        return a == b
    return type(a) is type(b) and a == b


def check(root, schema, value, path, errors):
    unknown = set(schema) - KNOWN
    if unknown:
        raise SchemaError(f"{path}: unsupported keyword(s) {sorted(unknown)}")
    if "$ref" in schema:
        check(root, resolve(root, schema["$ref"]), value, path, errors)
    for sub in schema.get("allOf", []):
        check(root, sub, value, path, errors)
    if "type" in schema:
        types = schema["type"] if isinstance(schema["type"], list) else [schema["type"]]
        bad = [t for t in types if t not in TYPES]
        if bad:
            raise SchemaError(f"{path}: unknown type {bad[0]}")
        if not any(TYPES[t](value) for t in types):
            errors.append(f"{path}: expected {'/'.join(types)}, got {type(value).__name__} {value!r:.60}")
            return
    if "enum" in schema and not any(same(value, e) for e in schema["enum"]):
        errors.append(f"{path}: {value!r:.60} is not one of {schema['enum']}")
    if "minimum" in schema and TYPES["number"](value):
        if isinstance(value, float) and math.isnan(value):
            errors.append(f"{path}: nan is not a number >= minimum {schema['minimum']}")
        elif value < schema["minimum"]:
            errors.append(f"{path}: {value} < minimum {schema['minimum']}")
    if isinstance(value, dict):
        for key in schema.get("required", []):
            if key not in value:
                errors.append(f"{path}: missing required key '{key}'")
        props = schema.get("properties", {})
        for key, sub in props.items():
            if key in value:
                check(root, sub, value[key], f"{path}.{key}", errors)
        extra = schema.get("additionalProperties", True)
        for key in value:
            if key in props:
                continue
            if extra is False:
                errors.append(f"{path}: key '{key}' is not allowed (additionalProperties: false)")
            elif isinstance(extra, dict):
                check(root, extra, value[key], f"{path}.{key}", errors)
    if isinstance(value, list):
        if "minItems" in schema and len(value) < schema["minItems"]:
            errors.append(f"{path}: {len(value)} items < minItems {schema['minItems']}")
        if "maxItems" in schema and len(value) > schema["maxItems"]:
            errors.append(f"{path}: {len(value)} items > maxItems {schema['maxItems']}")
        if "items" in schema:
            for i, item in enumerate(value):
                check(root, schema["items"], item, f"{path}[{i}]", errors)


def validate(schema, doc):
    """The errors of doc against schema, in og dev schema's order and words."""
    errors = []
    check(schema, schema, doc, "$", errors)
    return errors
