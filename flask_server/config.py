import os


class Config:
    SECRET_KEY = 'secret key'
    SQLALCHEMY_COMMIT_ON_TEARDOWN = True

    @staticmethod
    def init_app(app):
        pass


class DevelopmentConfig(Config):
    debug = True
    SQLALCHEMY_DATABASE_URI = os.environ.get('GENOME_EDITING_URI')



class TestingConfig(Config):
    pass


class ProductionConfig(Config):
    pass


config = {
    'development': DevelopmentConfig,
    'testing': TestingConfig,
    'production': ProductionConfig,
    'default': DevelopmentConfig
}
