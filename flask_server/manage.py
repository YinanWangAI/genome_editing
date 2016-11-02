#! usr/bin/env python
import os
from app import create_app, db
from flask_script import Manager, Server
from flask_migrate import Migrate, MigrateCommand
from flask_script import Shell
from app.models import SgrnaDesign

# import sys
# sys.path.append('/Users/yinan/PycharmProjects/')


app = create_app(os.getenv('FLASK_CONFIG') or 'default')
app.debug = True
manager = Manager(app)
migrate = Migrate(app, db)
manager.add_command('db', MigrateCommand)
manager.add_command("runserver", Server())
server = Server(host="0.0.0.0", port=8080)


def make_shell_context():
    return dict(app=app, db=db, SgrnaDesign=SgrnaDesign)
manager.add_command("shell", Shell(make_context=make_shell_context))

if __name__ == '__main__':
    manager.run()
