import json, sys
fn = sys.argv[1]
fam = sys.argv[2] if len(sys.argv) > 2 else None
J = json.load(open(fn))


def show(x, depth=0, maxd=4, key=''):
    pad = '  ' * depth
    if isinstance(x, dict):
        print(f'{pad}{key} dict[{len(x)}] keys={list(x.keys())[:25]}')
        if depth < maxd:
            for k in list(x.keys())[:25]:
                v = x[k]
                if isinstance(v, (dict, list)):
                    show(v, depth + 1, maxd, k)
                else:
                    print(f'{pad}  {k} = {str(v)[:150]}')
    elif isinstance(x, list):
        print(f'{pad}{key} list[{len(x)}]')
        if x and depth < maxd:
            show(x[0], depth + 1, maxd, key + '[0]')
    else:
        print(f'{pad}{key} = {str(x)[:150]}')


if isinstance(J, dict):
    print('top keys', list(J.keys())[:40])
    if fam and fam in J:
        show(J[fam], 0, 3, fam)
    elif fam:
        for k, v in J.items():
            if fam in str(k):
                show(v, 0, 3, k)
                break
    else:
        k0 = list(J.keys())[0]
        show(J[k0], 0, 3, k0)
else:
    print('top list', len(J))
    for item in J:
        if fam is None or (isinstance(item, dict) and fam in str(item.get('name', item.get('label', '')))):
            show(item, 0, 3, 'item')
            break
