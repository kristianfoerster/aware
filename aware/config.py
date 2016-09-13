from __future__ import print_function, division, absolute_import
import munch
import pprint

class Config(munch.Munch):
    @classmethod
    def from_pyfile(self, filename):
        conf = self()

        with open(filename) as f:
            tmp = {}
            exec(compile(f.read(), filename, 'exec'), tmp)
            conf.update({k: v for k, v in tmp.items() if k != '__builtins__'})

        return conf

    def __repr__(self):
        return pprint.pformat(dict(self))
