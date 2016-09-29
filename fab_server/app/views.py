import sys
sys.path = ['/Users/yinan/PycharmProjects/'] + sys.path

from flask import render_template, request, url_for, redirect, send_from_directory
from flask.ext.appbuilder import AppBuilder, BaseView, expose, has_access
from . import appbuilder

from flask.ext.appbuilder import ModelView
from flask.ext.appbuilder.models.sqla.interface import SQLAInterface

import genome_editing.design_sgRNA.design


class ToolsView(BaseView):

    default_view = 'design_sgrna'
    route_base = "/Tools"

    @expose('/Design sgRNA/', methods=['GET'])
    # @has_access
    def design_sgrna(self):
        # do something with param1
        # and return to previous page or index
        self.render_template('design_sgrna_main.html',
                             main_title='Design sgRNA')
        return self.render_template('design_sgrna_main.html',
                                    main_title='Design sgRNA')

    @expose('/Design sgRNA/', methods=['POST'])
    def catch_input(self):
        self.gene_info = request.form['gene_info']
        self.gene_info = [x.strip() for x in self.gene_info.split('\n')]
        self.pam = request.form['PAM']
        self.upstream_length = int(request.form['upstream_length'])
        self.downstream_length = int(request.form['downstream_length'])
        self.filter_tttt = request.form['filter_tttt']
        self.sgrna_length = int(request.form['sgrna_length'])
        return redirect(url_for('ToolsView.design_sgrna_result'))

    @expose('/Design sgRNA/result')
    def design_sgrna_result(self):
        if self.filter_tttt == 'Yes':
            filter_tttt = True
        else:
            filter_tttt = False

        flag = True
        for entrez_id in self.gene_info:
            entrez_id = int(entrez_id)
            sgrna_designer = genome_editing.design_sgRNA.design.Designer(
                entrez_id=entrez_id, sgrna_upstream=self.upstream_length,
                sgrna_downstream=self.downstream_length,
                sgrna_length=self.sgrna_length, filter_tttt=filter_tttt
            )
            sgrna_designer.get_sgrnas(pams=[self.pam])
            if flag:
                sgrna_designer_out = sgrna_designer.output()
                flag = False
            else:
                sgrna_designer_out = sgrna_designer_out.append(
                    sgrna_designer.output())
        sgrna_designer_out.to_csv('./results/design_sgrna_output.csv',
                                  index=None)
        return self.render_template('design_sgrna_result.html')

    @expose('/Score sgRNA/')
    # @has_access
    def score_sgrna(self):
        # do something with param1
        # and render template with param
        return self.render_template('score_sgrna_main.html',
                                    main_title='Under Development...')

    @expose('/Generate Negative Controls/')
    # @has_access
    def generate_negative_controls(self):
        # do something with param1
        # and render template with param
        return self.render_template('generate_nc_main.html',
                                    main_title='Under Development...')


appbuilder.add_view(ToolsView, "Design sgRNA", category='Tools')
appbuilder.add_link("Score sgRNA", href='/Tools/Score sgRNA/',
                    category='Tools')
appbuilder.add_link("Generate Negative Controls",
                    href='/Tools/Generate Negative Controls/',
                    category='Tools')


# class GroupModelView(ModelView):
#     datamodel = SQLAInterface(ContactGroup)
#     related_views = [ContactModelView]
